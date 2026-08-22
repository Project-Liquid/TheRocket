import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import CoolProp.CoolProp as CP

class Tank:
        
    @staticmethod
    def psi_to_Pa(psi):
        return psi * 6894.76

    @staticmethod
    def Pa_to_psi(Pa):
        return Pa / 6894.76

    @staticmethod
    def f_to_k(f):
        return (f - 32) * 5/9 + 273.15

    @staticmethod
    def k_to_f(k):
        return (k - 273.15) * 9/5 + 32

    def __init__(
            self, outlet_cv, inlet_cv, 
            outer_diameter, wall_thickness, max_height, outlet_height, 
            air_convection, fluid_convection, thermal_conductivity, 
            outlet_pressure, inlet_pressure, ambient_temp
            ):
        self.outlet_cv = outlet_cv
        self.inlet_cv = inlet_cv
        self.outer_diameter = outer_diameter
        self.wall_thickness = wall_thickness
        self.max_height = max_height
        self.outlet_height = outlet_height
        self.air_convection = air_convection
        self.fluid_convection = fluid_convection
        self.thermal_conductivity = thermal_conductivity
        self.outlet_pressure = outlet_pressure
        self.inlet_pressure = inlet_pressure
        self.ambient_temp = ambient_temp
        self.temp = self.ambient_temp
        self.fill_quality = 0
        self.drain_quality = 1
        self.name = "Tank"
        self.feed_tank_attached = False
        self.mass_flow_fill = 0
        self.mass_flow_drain = 0
        self.smoothing_factor = 0.3  # smoothing factor, lower = more damping
    
        self.cross_section_area = np.pi * (self.outer_diameter/2 - self.wall_thickness)**2
        self.tank_height = self.max_height
        self.tank_volume = self.cross_section_area * self.tank_height

        # Calculate heat transfer constants
        outer_surface_area = 2 * np.pi * self.outer_diameter/2 * self.tank_height
        inner_surface_area = 2 * np.pi * (self.outer_diameter/2 - self.wall_thickness) * self.tank_height
        mean_surface_area = (outer_surface_area - inner_surface_area) / np.log(outer_surface_area / inner_surface_area)
        self.thermal_conductance = 1/(1/(self.air_convection*outer_surface_area) + self.wall_thickness/(self.thermal_conductivity*mean_surface_area) + 1/(self.fluid_convection*inner_surface_area))

        # Convergence parameters
        self.energy_threshold = 1e-4 # Maximum energy differential at temperature convergence
        self.mass_threshold = 0.01 # kg # Maximum mass for small time steps near convergence
        self.min_dp = 0 # Pa # Pressure differential deadzone 
        self.converged = False
        self.in_error_state = False

        # Prepare data logging
        self.data = []
        self._prev_mass_flow_fill = 0
        self._prev_mass_flow_drain = 0


    def fetch_properties(self):
        self.liquid_density = CP.PropsSI('D', 'T', self.temp, 'Q', 0, 'NitrousOxide')
        self.vapor_density = CP.PropsSI('D', 'T', self.temp, 'Q', 1, 'NitrousOxide')
        self.liquid_energy = CP.PropsSI('U', 'T', self.temp, 'Q', 0, 'NitrousOxide')
        self.vapor_energy = CP.PropsSI('U', 'T', self.temp, 'Q', 1, 'NitrousOxide')
        self.tank_pressure = CP.PropsSI('P', 'T', self.temp, 'Q', 0, 'NitrousOxide')

    def set_fill_height(self, temp, liquid_height):
        self.temp = temp
        self.liquid_height = liquid_height
        self.vapor_height = self.tank_height - self.liquid_height
        self.fetch_properties()
        self.liquid_mass = self.liquid_height * self.cross_section_area * self.liquid_density
        self.vapor_mass = self.vapor_height * self.cross_section_area * self.vapor_density
        self.total_mass = self.liquid_mass + self.vapor_mass
        self.energy = self.liquid_mass * self.liquid_energy + self.vapor_mass * self.vapor_energy

    def set_fill_mass(self, temp, total_mass):
        self.temp = temp
        self.fetch_properties()
        self.total_mass = total_mass
        # Calculate liquid and vapor masses based on total mass and densities
        self.liquid_mass = (self.tank_volume - self.total_mass / self.vapor_density) * (self.liquid_density * self.vapor_density / (self.vapor_density - self.liquid_density))
        self.vapor_mass = self.total_mass - self.liquid_mass
        # Update heights based on masses and densities
        self.liquid_height = self.liquid_mass / (self.cross_section_area * self.liquid_density)
        self.vapor_height = self.tank_height - self.liquid_height
        # Update energy based on new masses and energies
        self.energy = self.liquid_mass * self.liquid_energy + self.vapor_mass * self.vapor_energy
  

    def update_mass(self, time_step):
        self.liquid_height = self.liquid_mass / (self.cross_section_area * self.liquid_density)
        self.drain_quality = 0 if self.liquid_height > self.outlet_height else 1
        fill_density = self.vapor_density * self.fill_quality + self.liquid_density * (1 - self.fill_quality)
        drain_density = self.vapor_density * self.drain_quality + self.liquid_density * (1 - self.drain_quality)
        if drain_density == 0:
            drain_density = 1e-6  # Avoid division by zero
        
        if not self.feed_tank_attached:
            self.mass_flow_fill = fill_density * self.inlet_cv * abs(self.Pa_to_psi(self.inlet_pressure - self.tank_pressure) * 1000 / fill_density) ** 0.5 * 0.00379 / 60
        self.mass_flow_drain = drain_density * self.outlet_cv * abs(self.Pa_to_psi(self.tank_pressure - self.outlet_pressure) * 1000 / drain_density) ** 0.5 * 0.00379 / 60
        # Reverse flow if pressure gradient flips
        if self.inlet_pressure < self.tank_pressure: 
            self.mass_flow_fill *= -1
        if self.tank_pressure < self.outlet_pressure:
            self.mass_flow_drain *= -1

        # Apply deadzone for small pressure differentials
        if abs(self.inlet_pressure - self.tank_pressure) < self.min_dp:
            self.mass_flow_fill = 0
        if abs(self.tank_pressure - self.outlet_pressure) < self.min_dp:
            self.mass_flow_drain = 0

        # Smooth mass flow rates to avoid oscillations
        if not self.feed_tank_attached:
            self.mass_flow_fill = self.smoothing_factor * self.mass_flow_fill + (1 - self.smoothing_factor) * self._prev_mass_flow_fill
        self.mass_flow_drain = self.smoothing_factor * self.mass_flow_drain + (1 - self.smoothing_factor) * self._prev_mass_flow_drain

        mass_flow = self.mass_flow_fill - self.mass_flow_drain
        runout_adjustment = 1.0
        if -mass_flow * time_step > self.total_mass:
            runout_adjustment = abs(self.total_mass / (mass_flow * time_step))
        dm = mass_flow * time_step * runout_adjustment
        self.total_mass += dm
        self.dm_liquid_in = self.mass_flow_fill * (1-self.fill_quality) * time_step * runout_adjustment
        self.dm_vapor_in = self.mass_flow_fill * self.fill_quality * time_step * runout_adjustment
        self.dm_liquid_out = -self.mass_flow_drain * (1-self.drain_quality) * time_step * runout_adjustment
        self.dm_vapor_out = -self.mass_flow_drain * self.drain_quality * time_step * runout_adjustment
        # print(f"{self.name}: dp={self.Pa_to_psi(self.tank_pressure - self.outlet_pressure)}, runout_adjustment={runout_adjustment}")

    def update_temperature(self, time_step, time):
        # Calculate change in energy
        target_energy = self.energy
        heat_transfer = self.thermal_conductance * (self.ambient_temp - self.temp) * time_step
        target_energy += heat_transfer
        #prev_temp = self.temp

        # Liquid and vapor
        if (self.liquid_mass > 0 or self.mass_flow_fill > 0):
            h_liquid = CP.PropsSI('H', 'T', self.temp, 'Q', 0, 'NitrousOxide')
            h_vapor = CP.PropsSI('H', 'T', self.temp, 'Q', 1, 'NitrousOxide')

            # Inlet fluid enthalpy at supply pressure (saturated liquid)
            if self.inlet_cv > 0 and self.inlet_pressure > 0 and self.inlet_pressure < 7.244e6:
                h_inlet = CP.PropsSI('H', 'P', self.inlet_pressure, 'Q', 0, 'NitrousOxide')
            else:
                h_inlet = 0

            # Energy balance: incoming fluid + outgoing fluid
            inlet_energy = self.dm_liquid_in * h_inlet if self.mass_flow_fill > 0 else 0
            outlet_energy = self.dm_liquid_out * h_liquid + self.dm_vapor_out * h_vapor
            target_energy += inlet_energy + outlet_energy
        
            # Find temperature and corresponding vapor quality with target energy
            last_temp = self.temp
            last_energy = self.energy
            for i in range(1000):
                error = target_energy - self.energy
                if abs(error) <= self.energy_threshold:
                    break

                # Newton's method
                # Estimate dU_dT
                dT = 0.01
                rho_l_plus = CP.PropsSI('D', 'T', self.temp + dT, 'Q', 0, 'NitrousOxide')
                rho_v_plus = CP.PropsSI('D', 'T', self.temp + dT, 'Q', 1, 'NitrousOxide')
                drho_l_dT = (rho_l_plus - self.liquid_density) / dT
                drho_v_dT = (rho_v_plus - self.vapor_density) / dT
                # m_l = (V - m/rho_v) * rho_l*rho_v / (rho_v - rho_l)
                # Let A = V - m/rho_v, B = rho_l*rho_v/(rho_v - rho_l)
                A = self.tank_volume - self.total_mass / self.vapor_density
                B = self.liquid_density * self.vapor_density / (self.vapor_density - self.liquid_density)
                
                dA_dT = self.total_mass / self.vapor_density**2 * drho_v_dT
                dB_dT = ((drho_l_dT * self.vapor_density + self.liquid_density * drho_v_dT) * (self.vapor_density - self.liquid_density)
                        - self.liquid_density * self.vapor_density * (drho_v_dT - drho_l_dT)
                    ) / (self.vapor_density - self.liquid_density)**2
                
                dm_l_dT = dA_dT * B + A * dB_dT
                dm_v_dT = -dm_l_dT
                u_l = CP.PropsSI('U', 'T', self.temp, 'Q', 0, 'NitrousOxide')
                u_v = CP.PropsSI('U', 'T', self.temp, 'Q', 1, 'NitrousOxide')
                cv_l = CP.PropsSI('CVMASS', 'T', self.temp, 'Q', 0, 'NitrousOxide')
                cv_v = CP.PropsSI('CVMASS', 'T', self.temp, 'Q', 1, 'NitrousOxide')
                
                dU_dT = (self.liquid_mass * cv_l 
                    + self.vapor_mass * cv_v 
                    + (u_v - u_l) * dm_v_dT)
                #slope = -5e2
                slope = dU_dT
                if (self.energy != last_energy and self.temp != last_temp):
                    slope = (self.energy - last_energy) / (self.temp - last_temp)
                last_energy = self.energy
                last_temp = self.temp
                self.temp += error / slope
                self.temp = min(self.temp, 309.3)
                self.temp = max(self.temp, 182.4)
                
                self.liquid_density = CP.PropsSI('D', 'T', self.temp, 'Q', 0, 'NitrousOxide')
                self.vapor_density = CP.PropsSI('D', 'T', self.temp, 'Q', 1, 'NitrousOxide')
                self.liquid_energy = CP.PropsSI('U', 'T', self.temp, 'Q', 0, 'NitrousOxide')
                self.vapor_energy = CP.PropsSI('U', 'T', self.temp, 'Q', 1, 'NitrousOxide')
                
                self.liquid_mass = (self.tank_volume - self.total_mass / self.vapor_density) * (self.liquid_density * self.vapor_density / (self.vapor_density - self.liquid_density))
                self.liquid_mass = min(max(self.liquid_mass, 0), self.total_mass) # Clamp to [0, total_mass]
                self.vapor_mass = self.total_mass - self.liquid_mass
                self.energy = self.liquid_mass * self.liquid_energy + self.vapor_mass * self.vapor_energy
                
            else:
                print("Warning: energy convergence failed at t=%.2f" % time)
                print("Target: ", target_energy)
                print("energy: ", self.energy)
                print("temp: ", self.temp)
                
            self.energy = target_energy
            # Update pressure at new temperature
            self.tank_pressure = CP.PropsSI('P', 'T', self.temp, 'Q', 1, 'NitrousOxide')
            # Transition to vapor-only
            if self.liquid_mass <= 0:
                print("Liquid evaporated in %.1fs" % time)
                self.vapor_mass = self.total_mass
                self.liquid_mass = 0
                self.vapor_density = self.total_mass / self.tank_volume
                u = self.energy / self.total_mass
                
                # At evaporation point, energy is already correctly calculated as internal energy
                # Just recalculate T and P for the pure vapor state
                try:
                    self.temp = CP.PropsSI('T', 'D', self.vapor_density, 'U', u, 'NitrousOxide')
                    self.tank_pressure = CP.PropsSI('P', 'D', self.vapor_density, 'U', u, 'NitrousOxide')
                except ValueError:
                    print(f"Warning: Invalid state at evaporation point")
                    # Keep previous values and let vapor-only section handle it
                    pass
                
        # Vapor only 
        else:
            self.vapor_mass = self.total_mass
            self.liquid_mass = 0
            
            # Simple energy balance for open system (tank losing mass):
            # dE/dt = Q - h_out * dm_out/dt
            # where h_out is specific enthalpy at tank conditions
            
            # Apply heat transfer and mass flow
            vapor_density_old = self.total_mass / self.tank_volume
            
            try:
                # Get enthalpy at current state to subtract from exiting mass
                h_vapor_current = CP.PropsSI('H', 'T', self.temp, 'D', vapor_density_old, 'NitrousOxide')
                
                # Heat transfer
                heat_transfer = self.thermal_conductance * (self.ambient_temp - self.temp) * time_step

                # Energy balance: energy change = heat in + enthalpy of exiting mass
                # (dm_vapor_out is negative when mass exits, so we subtract it as in dE = Q - h*|dm_out|)
                self.energy = self.energy + heat_transfer + self.dm_vapor_out * h_vapor_current
              
                # Now find temperature at the new energy state
                u_specific = self.energy / self.total_mass
                vapor_density_new = self.total_mass / self.tank_volume
                
                # Solve for temperature using internal energy
                self.temp = CP.PropsSI('T', 'D', vapor_density_new, 'U', u_specific, 'NitrousOxide')
                self.tank_pressure = CP.PropsSI('P', 'D', vapor_density_new, 'U', u_specific, 'NitrousOxide')
                
            except ValueError as e:
                print(f"Warning: Invalid vapor state at t={time:.2f}, terminating")
                print(f"  Reason: {str(e)}")
                self.in_error_state = True

    def get_pressure_psi(self):
        return self.Pa_to_psi(self.tank_pressure)
    
    def get_temp_f(self):
        return self.k_to_f(self.temp)
    
    def feed_from_tank(self, other_tank):
        self.feed_tank_attached = True
        self.inlet_pressure = other_tank.tank_pressure
        self.fill_quality = other_tank.drain_quality
        self.mass_flow_fill = other_tank.mass_flow_drain

    def drain_to_tank(self, other_tank):
        self.outlet_pressure = other_tank.tank_pressure

    def log_state(self, time):
        self.data.append({
            "Time": time,
            "Pressure": self.tank_pressure,
            "Temperature": self.temp,
            "Mass": self.total_mass,
            "Liquid": self.liquid_mass,
            "Energy": self.energy,
            "Vapor Flow In": self.dm_vapor_in,
            "Liquid Flow In": self.dm_liquid_in,
            "Vapor Flow Out": self.dm_vapor_out,
            "Liquid Flow Out": self.dm_liquid_out,
        })

    def save_state(self):
        self._prev_pressure    = self.tank_pressure
        self._prev_temp        = self.temp
        self._prev_mass        = self.total_mass
        self._prev_mass_flow_fill = self.mass_flow_fill
        self._prev_mass_flow_drain = self.mass_flow_drain
        self._prev_liquid_mass = self.liquid_mass
        self._prev_energy      = self.energy

