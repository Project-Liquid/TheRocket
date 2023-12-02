from geopy import distance, Point
from googlemaps import addressvalidation, client, directions, distance_matrix, elevation, exceptions, geocoding, geolocation, maps, places, roads, timezone
from math import log
import math
import random
import sys

class HFD:
    """Provides methods for calculating the HFD of a worst-case fuel and oxidizer rupture event."""
    def __init__(self) -> None:
        pass

    def __get_area(self, length, diameter, presented=True):
        """Returns presented area A_p"""
        area = 2 * math.pi * (diameter)/2 * (length) + (2 * math.pi * ((diameter)/2)**2)
        return 0.25 * area if presented else area
    



    def __get_m_0(self, mean_shape_factor, diameter, length, material_density, A_p=None):
        """Returns the mass m_0"""
        presented_area = self.__get_area(length, diameter) if A_p == None else A_p
        m_0 = (((presented_area/mean_shape_factor)**(1.5)) * material_density)
        return m_0




    def __get_V(self, M, C, A_p, V_g = 2.44):
        """Returns the gurney velocity V in mm per µs."""
        return (V_g)/(((M/C)+(A_p))**0.5)
    



    def __get_v(self, V, R, L):
        """Returns the velocity v at a given distance along the trajectory."""
        return V * math.exp(-R / L)
    



    def __get_L1(self, k, C_d=1.28, rho=1.225):
        """Returns the value for length L_1.
        Args:
            k: Ballastic density in g/cm^3
            C_d: Drag coefficient
            rho: The atmopsheric density
        """
        #x100 because we convert from g/cm^3 to gcm, then multiply by 100 to get meters
        return (2 * (((k)**2)**(1/3)))/(C_d * rho)




    def __get_L(self, m, L_1):
        """Return the value of length L, where the distance fragment velocity is 1/e its value."""
        return L_1 * (m**(1/3))
    



    def __get_m(self, R, L_1, E_cr, V, m_last = None, m_ter = None, m_raise=True, m_iter=0.01):
        """Returns the value for m, the minimum mass that shall be classified as a hazardous fragement. The smallest of two values, m1 and m2, are returned."""
        # Consider m_2
        m_2 = ((2.0 * E_cr) / (9.81 * L_1)) ** (3/4) if m_ter == None else m_ter
        A = 2*E_cr
        q = lambda m : m*(V*1000**2)*math.exp((-2.0*R)/(L_1*(m**(1.0/3.0))))

        # for i in range(1, 989):
        #     m = (0.9 + (i * 0.1))*0.001
        #     a_iter = a(m)
        #     diff = A - a_iter
        #     if abs(diff) < 0.5:
        #         return min(m, m_2)        
        # return None

        m_new = m_last if m_last != None else m_2
        m_new += m_iter if m_raise else -m_iter

        a_new = q(m_new)
        a_2 = q(m_2/100)
        a_diff = a_new - A
        try:
            if (abs(a_diff) < 0.5):
                print(m_new)
                print(m_2)
                return min(m_new*100, m_2)
            elif (abs(a_2 - A) < 0.5):
                print("asdfasdfasdfasdfadsf")
                return min(m_new*100, m_2)

            a_old = q(m_last if m_last != None else m_2)

            if round(a_new, 5) == round(a_old, 5):
                return None
            # elif a_new < 0.000001:
            #     return None
            
            condit = (abs(a_new) > abs(a_old)) or (a_diff <= 0 and a_new > 0)
            return self.__get_m(R, L_1, E_cr, V, m_last=m_new, m_ter=m_ter, m_raise=abs(m_raise - condit), m_iter=m_iter if not condit else m_iter*0.1)
        except Exception as err:
            if type(err) == TypeError:
                return self.__get_m(R, L_1, E_cr, V, m_last=m_new, m_ter=m_ter, m_raise=abs(m_raise - 1), m_iter=m_iter*0.1)
            return None




    def __get_N(self, M_t, m_0, m):
        """Returns the number of fragements N."""
        return (M_t/m_0)*math.exp(-((2*m)/m_0)**0.5)




    def __get_Q_0(self, N, divisions=4.0):
        """Returns the number of fragements per solid unit angle."""
        return N/divisions
    



    def __get_c_as_tnt(self, C, key, grams_per_mole=None, kj_per_mole=None):
        """For a given mass of chemical C in kg, returns the equivalent TNT mass in kg. Does not consider yield."""
        chemicals = {
            "ethane": [30.07, 1560.0],
            "nitrous": [44.01, 82.0]
        }

        is_valid = key in chemicals.keys()

        if (not is_valid):
            assert(grams_per_mole != None and kj_per_mole != None)

        grams = C * 1000.0
        c_grams_per_mole = grams / (chemicals[key][0] if is_valid else grams_per_mole)
        c_kj_per_mole = c_grams_per_mole * (chemicals[key][1] if is_valid else kj_per_mole)
        kg_tnt = c_kj_per_mole / (4.184*(10**3))
        return kg_tnt




    def __get_r_and_m(self, L_1, E_cr, V, m_0, M_t, R, target_p=0.01, target_q=(1/55.7418), r_iter=1, A_t=0.58, divisions=4.0):
        """Considers the values of r and m for a given R"""
        m = self.__get_m(R, L_1, E_cr, V)
        if m == None:
            return None
            # m = ((2.0 * E_cr) / (9.81 * L_1)) ** (3/4)
            # a = lambda m : m*1000*(V**2)*math.exp((-2.0*R)/(L_1*(m*1000**(1.0/3.0))))
            # print(a(m=m))

        N = self.__get_N(M_t, m_0, m)
        q_0 = N/4
        q = (q_0/(R**2.0))*math.exp(-(2.0*m/m_0)**0.5)
        l = self.__get_L(m, L_1)
        
        if ((1 - math.exp(-q * A_t)) <= target_p) and (q <= target_q):
            print("Q: " + str(q))
            return R, m
        else:
            return None

    def get_radius_and_m(self, ethane_tank_mass, ethane_mass, nitrous_tank_mass, nitrous_mass, phi=None, vg = 2.44, c_d = 1.28, rho=1.225, E_cr=79.0, A=0.5, target_p=0.01, target_q=(1.0/55.7418), r_start=1.0, r_iter=1.0, A_t=0.58, divisions=4.0, yield_factor=1.0, fuel_tank_diameter=0.18, fuel_tank_length=0.84, fuel_tank_density=2.710, oxidizer_tank_diameter=0.18, oxidizer_tank_length=0.41, oxidizer_tank_density=2.710, only_tanks=True):
        """Returns the radius (meters) and the minimum fragement mass (kg).
        Args:
            ethane_tank_mass: The dry mass of the ethane tank in kg.
            ethane_mass: The mass of the ethane within the tank.
            nitrous_tank_mass: The dry mass of the nitrous tank in kg.
            nitrous_mass: The mass of the nitrous within the tank.
            k: The ballistic density and/or shape factor of the average fragement (g/cm^3). It is recommended to be either 2.37 or conservatively 4.74.
            vg: The gurney constant for the equivalent explosive. Here, the propellants are equivalated for TNT equivalency, which is 2.44.
            c_d: The drag coefficient of the average fragement.
            rho: The atmospheric density
            E_cr: The maximum tolerable energy of a hazardous fragement at a given distance. The DoD has stated that 79 Joules or 58lb-ft is the limit.
            A: 0.5, derived from n/n+2, where 2 is used for a cylindrical shape.
            target_p: The probability value of p for which 1-exp(-q*A_t) is less than or equal too.
            target_q: The areal density. The DoD suggests this area is 1/600ft^2, or 1/55.7418m^2
            r_start: The first radius value (meters) to start calculating. Default is 1m.
            r_iter: The value in meters by which to check the radius r. Default is 1m.
            A_t: The target area A_t, for someone facing the explosion and not taking cover. The default is recommended as 0.58m^2 by technical document 12, although this is a very generous use case.
            divisions: The standard number of solid angle divisions used to calculate Q_0.
            yield_factor: The yield of the explosive mass, if applicable between 0.0 and 1.0. Default is 1.0.
        Returns:
            radius: The radius of the HFD in meters
            m: The minimum fragement mass (kg).
        """
        chem_dict = {
            "ethane": [ethane_tank_mass, ethane_mass],
            "nitrous": [nitrous_tank_mass, nitrous_mass]
        }

        # if rocket_sim:
        #     k_tanks, sf_tanks = self.get_ballistic_density_from_sf([(sf_fuel+sf_oxidizer)/2], material_density=max(fuel_tank_density, oxidizer_tank_density), return_shape_factors=True)#min(k_fuel + k_oxidizer), max(sf_fuel+sf_oxidizer)
        #     k_tanks, sf_tanks = k_tanks[0], sf_tanks[0]
        #     R_tanks, m_tanks = self.get_r_and_m_from_chem_dict(chem_dict, k_tanks, vg, c_d, rho, E_cr, A, target_p, target_q, r_start+r_iter, r_iter, A_t, divisions, yield_factor, fuel_tank_length + oxidizer_tank_length, (fuel_tank_diameter + oxidizer_tank_diameter)/2, max(fuel_tank_density, oxidizer_tank_density), sf_tanks)
        #     return R_tanks, m_tanks

        # Get factor for fuel
        k_fuel, sf_fuel = self.get_ballistic_density(fuel_tank_length, fuel_tank_diameter, values=[phi] if phi != None else ["mean"], material_density=fuel_tank_density, return_shape_factors=True)
        k_fuel, sf_fuel = k_fuel[0], sf_fuel[0]
        if not only_tanks:
            R_fuel, m_fuel = self.get_r_and_m_from_chem_dict({"ethane": [ethane_tank_mass, ethane_mass]}, k_fuel, vg, c_d, rho, E_cr, A, target_p, target_q, r_start+r_iter, r_iter, A_t, divisions, yield_factor, fuel_tank_length, fuel_tank_diameter, fuel_tank_density, sf_fuel)
            print("[Fuel]: Radius=" + str(R_fuel) + ";m=" + str(m_fuel))

        # Get oxidizer
        k_oxidizer, sf_oxidizer = self.get_ballistic_density(oxidizer_tank_length, oxidizer_tank_diameter, values=[phi] if phi != None else ["mean"], material_density=oxidizer_tank_density, return_shape_factors=True)
        k_oxidizer, sf_oxidizer = k_oxidizer[0], sf_oxidizer[0]
        if not only_tanks:
            R_oxidizer, m_oxidizer = self.get_r_and_m_from_chem_dict({"nitrous": [nitrous_tank_mass, nitrous_mass]}, k_oxidizer, vg, c_d, rho, E_cr, A, target_p, target_q, r_start+r_iter, r_iter, A_t, divisions, yield_factor, oxidizer_tank_length, oxidizer_tank_diameter, oxidizer_tank_density, sf_oxidizer)
            print("[Ox]: Radius=" + str(R_oxidizer) + ";m=" + str(m_oxidizer))


        # Get tanks
        #k_tanks, sf_tanks = self.get_ballistic_density_from_sf([(sf_fuel+sf_oxidizer)/2], material_density=max(fuel_tank_density, oxidizer_tank_density), return_shape_factors=True)#min(k_fuel + k_oxidizer), max(sf_fuel+sf_oxidizer)
        k_tanks, sf_tanks =  self.get_ballistic_density_from_sf([min(sf_fuel, sf_oxidizer)], material_density=max(fuel_tank_density, oxidizer_tank_density), return_shape_factors=True)
        k_tanks, sf_tanks = k_tanks[0], sf_tanks[0]
        R_tanks, m_tanks = self.get_r_and_m_from_chem_dict(chem_dict, k_tanks, vg, c_d, rho, E_cr, A, target_p, target_q, r_start+r_iter, r_iter, A_t, divisions, yield_factor, fuel_tank_length + oxidizer_tank_length, (fuel_tank_diameter + oxidizer_tank_diameter)/2, max(fuel_tank_density, oxidizer_tank_density), sf_tanks)
        print("[All]: Radius=" + str(R_tanks) + ";m=" + str(m_tanks))

        return R_tanks, m_tanks
    



    def get_r_and_m_from_chem_dict(self, chem_dict, k=4.74, vg = 2.44, c_d = 1.28, rho=1.225, E_cr=79.0, A=0.5, target_p=0.01, target_q=(1/55.7418), r_start=1, r_iter=1, A_t=0.58, divisions=4.0, yield_factor=1.0, length=0.84, diameter=0.18, density=2.710, sf=1.384):
        """Returns the radius using the chem_dict"""
        m_total = 0
        c = 0
        divisions = 0
        for key, chem in chem_dict.items():
            temp_m = chem[0]
            temp_c = chem[1]
            tnt_c = self.__get_c_as_tnt(temp_c, key)
            m_total += temp_m
            c += tnt_c
            divisions += 1

        c = c * yield_factor

        A_p = self.__get_area(length, diameter)
        area = self.__get_area(length, diameter, False)
        v_gurney = self.__get_V(m_total, c, area, vg)
        l_1 = self.__get_L1(k, c_d, rho)
        m_0 = self.__get_m_0(sf, length=length, diameter=diameter, material_density=density, A_p=A_p)
        
        try:
            R_out, m_out = self.__get_r_and_m(l_1, E_cr, v_gurney, m_0, m_total, r_start, target_p, target_q, r_iter, A_t, divisions)
            return R_out, m_out
        except:
            if r_start == 600:
                print("Bruh")            
            return self.get_r_and_m_from_chem_dict(chem_dict, k, vg, c_d, rho, E_cr, A, target_p, target_q, r_start+r_iter, r_iter, A_t, divisions, yield_factor)
        
    def old_get_shape_factor(self, length, diameter):
        N = 1000  # default number of samples or specify on command line
        l_d = length/diameter  # default is a L/D = 1 cylinder, or specify on command line

        if len(sys.argv) == 2:
            N = int(sys.argv[1])  # number of samples
        elif len(sys.argv) == 3:
            l_d = float(sys.argv[1])  # L/D
            N = int(sys.argv[2])  # number of samples

        C = (math.pi / 4 * l_d) ** (-2 / 3)
        th, ph, sf = 0.0, 0.0, 0.0

        for i in range(N):
            th, ph = self.spherical_avoidance()
            sf = C * (l_d * math.sin(th) + math.pi / 4 * abs(math.cos(th)))


    def get_yaw_shape_factor(self, length, diameter, yaw):
        l_d = length/diameter
        sf = (math.pi/4.0 * (l_d))**(-2.0/3.0) * (l_d*math.sin(yaw) + math.pi/4.0 * abs(math.cos(yaw)))
        return sf
    
    def get_shape_factor(self, length, diameter, options):
        l_d = length/diameter
        return_options = []

        if "min" in options:
            min = self.get_yaw_shape_factor(length, diameter, 0 if l_d >= (math.pi/4.0) else math.pi/2.0)
            return_options.append(min)
        if "max" in options:
            max_yaw = math.atan(l_d/(math.pi/4.0))
            max = self.get_yaw_shape_factor(length, diameter, max_yaw)
            return_options.append(max)
        if "mean" in options:
            mean = (math.pi/4.0)**(1.0/3.0) * (l_d)**(-2.0/3.0) * (l_d + 0.5)
            return_options.append(mean)
        for option in options:
            if option not in ["min", "max", "mean"]:
                if type(option) == float or type(option) == int:
                    return_options.append(self.get_yaw_shape_factor(length, diameter, float(option) * math.pi / 180.0))

        return return_options
    
    def get_min_shape_factor(self, length, diameter):
        return self.get_shape_factor(length, diameter, ["min"])[0]
    
    def get_max_shape_factor(self, length, diameter):
        return self.get_shape_factor(length, diameter, ["max"])[0]
    
    def get_mean_shape_factor(self, length, diameter):
        return self.get_shape_factor(length, diameter, ["mean"])[0]
    
    def get_density_in_grains_per_cubic_inch(self, material_denstiy_lb_in3):
        #density * (453.592g/1lb) * (15.4324gr/1g)
        return material_denstiy_lb_in3 * 453.592 * 15.4324
    
    def get_ballistic_density(self, length, diameter, values, material_density=2.710, use_si_units=True, return_shape_factors=False):
        """
        Returns the ballistic density for expressed values of shape factor yaw, as both keywords and floats. Default returns in kg per m^3, otherwise in grains per inch^3.

        Args:
            length: length of the tank in units equal to diameter
            diameter: diameter of the tank in units equal to length
            material_density: density in grams per cm^3
            values: values to get density for. Can use string keywords "max", "min", and "mean", as well as any degree values for yaw between 0 and 180. 0 degrees is to be the front of the tank.
        """
        factors = self.get_shape_factor(length, diameter, values)
        k_factors = [252.9 * (material_density/(factor**(3/2))) * (61.0237/15.4324 if use_si_units else 1)  for factor in factors]
        if return_shape_factors:
            return k_factors, factors
        else:
            return k_factors
        
    def get_ballistic_density_from_sf(self, sf_factors, material_density, use_si_units=True, return_shape_factors=False):
        k_factors = [252.9 * (material_density/(factor**(3/2))) * (61.0237/15.4324 if use_si_units else 1)  for factor in sf_factors]
        if return_shape_factors:
            return k_factors, sf_factors
        else:
            return k_factors
        


class __BlastDistance:
    """NOT USED"""
    def __init__(self, fw_grams = 1300, ox_grams = 2000):
        self.fuel_weight = fw_grams / 453.6
        self.oxidizer_weight = ox_grams / 453.6

    def net_explosive_weight(self, input_fw = None, input_ow = None, yield_percentage = 0.2):
        return yield_percentage * ((self.fuel_weight if input_fw == None else input_fw) + (self.oxidizer_weight if input_ow == None else input_ow))
    
    def overpressure_distance(self, k_factor = 45, auto_w = True):
        w = self.net_explosive_weight() if auto_w == True else auto_w
        assert(type(w) == float or type(w) == int)
        return k_factor * (w**(1/3))
    
    def hfd_distance(self, auto_w = True):
        w = self.net_explosive_weight() if auto_w == True else auto_w
        assert(type(w) == float or type(w) == int)
        return -1133.9 + (389 * log(w))


class Blast:
    """Defines helpful methods for getting the blast radius"""
    def __init__(self, origin = Point(41.28803, -72.74378), key_client = "AIzaSyA3v9EoJh6RPagG491yLx-8cH9HYns2gKg"):
        """
        A class to determine to implement terrain-aware distance calculations.
        Args:
            origin: An origin of type Point(lat, long). Altitude should be none.
            key_client: The Google Maps client API key to use for calculations. If you use the default key, please be mindful of your usage (it costs about a bit more than a cent or so for an intensive run, but I get $200 in free credits per month).
        """
        self.origin = origin#Yes this is a lot of commented out points. Point(41.28804, -72.74374)#Point(41.288156, -72.743703)
        #Point(41.31298570804371, -72.92446707721605)#Point(41.288156, -72.743703) #Point(41.319463036328, -72.9203207745025)
        #origin = #Point(41.2878801, -72.7439151)#Point(41.2881304, -72.7437407)#Point(41.2879448, -72.7438069)#41°17'17"N 72°44'37"W
        self.key_client = client.Client(key_client)

    def get_distance(self, bearing, radius = 150, origin = None):
        new_coords = distance.distance(miles=(radius/5280.0)).destination(origin if origin != None else self.origin, bearing=bearing)
        return new_coords

    def get_distance_between_points(self, destination: Point, start: Point = False):
        return distance.geodesic(Point(self.origin[0], self.origin[1]).format_decimal(False) if start == False else start.format_decimal(False), destination.format_decimal(False)).ft
    
    def get_phi(self, phi, phi_offset=0):
        phi_prime = 180.0 - abs(phi-180+phi_offset)
        if phi_prime < 0:
            phi_prime += 360
        elif phi_prime >= 360:
            phi_prime -= 360
        return phi_prime


    def get_circular_points(self, radius = None, origin = None, num = 20, phi_offset = 0):
        frag = HFD()
        point_dict = {}
        if (360 % num != 0):
            print("Number is not divisible, will return empty points")
            return []
        
        for i in range(0, 360, int(360/num)):
            phi = self.get_phi(i, phi_offset)
            if radius == None:
                Radius, _ = frag.get_radius_and_m(
                    ethane_tank_mass=54.0,
                    ethane_mass=14.51,
                    nitrous_tank_mass=54.0,
                    nitrous_mass=27.2,
                    phi=abs(180.0 - abs(i-180+phi_offset)),
                    yield_factor=1.0,
                    fuel_tank_diameter=0.23,
                    fuel_tank_length=1.30,
                    fuel_tank_density=7.75,
                    oxidizer_tank_diameter=0.23,
                    oxidizer_tank_length=1.30,
                    oxidizer_tank_density=7.75,
                    only_tanks=True)
                Radius = Radius * 3.28084 # convert from meters to feet
            else:
                Radius = radius
            point_dict[i] = {
                "dist": self.get_distance(i, Radius, self.origin if origin == None else origin),
                "radius": Radius
            }
            print("Phi: " + str(abs(180.0 - abs(i-180+phi_offset))) + "; Radius: " + str(Radius))
        return point_dict

    def get_elevation_dict(self, point_dict, origin: Point, radius=None, division=25):
        path_points = {}
        for key, point in point_dict.items():
            elevated_points = elevation.elevation_along_path(client=self.key_client, path=[origin, point["dist"]], samples=division)
            chosen_point = self.get_best_elevated_point(elevated_points, point["radius"] if radius == None else radius)
            if chosen_point != None:
                path_points["" + str(point["dist"]) + "|" + str(self.get_distance_between_points(chosen_point, origin))] = chosen_point
        return path_points

    def get_point_from_elevated_point(self, elevated_point):
        return Point(latitude=elevated_point["location"]["lat"], longitude=elevated_point["location"]["lng"], altitude=(elevated_point["elevation"]/1000.0))


    def get_best_elevated_point(self, elevated_points, radius=150):
        origin = self.get_point_from_elevated_point(elevated_points[0])

        for i in range(1, len(elevated_points)):
            temp_point = self.get_point_from_elevated_point(elevated_points[i])
            ground_dist = self.get_distance_between_points(temp_point, origin)
            print(ground_dist)
            altitude = (temp_point.altitude - origin.altitude) * 3280.84
            print(altitude)
            ans = (((ground_dist**2)+((altitude)**2))**0.5) - radius
            if ans >= 0:
                print("Good at " + str(i))
                return temp_point
        
        print("Defaulting at " + str(len(elevated_points)))
        return temp_point

    def get_path_points_as_array(self, path_points):
        path_array = [value.format_decimal(False) for value in path_points.values()]
        path_array.append(path_array[0])
        return path_array

    def returnOffset(self, point, latFeet, longFeet):
        first_shift = distance.distance(miles=(latFeet/5280.0)).destination(point, bearing=0)
        return distance.distance(miles=(longFeet/5280.0)).destination(first_shift, bearing=90)

    def map_points(self, radius = None, input_origin = None, num = 30, division=25, latOffset = 0, longOffset = 0, is_map_centered_by_origin = False, label = "stony_creek", zoom=18, should_confirm=True, phi_offset=0):
        """Given a radius, return a map of the terrain-aware circular area. 
        Args:
            radius: The radius of the blast in feet.
            input_origin: The origin to use, otherwise, the origin defined in the class
            num: The number of points defining the perimiter of the circle. Default is 30.
            division: The number of data-mapped elevation points to consider for a given angle specified by num. A larger division yields more specific mapping information. Default is 25.
            latOffset: The north-south offset in feet. Default 0.
            longOffset: The west-east offset in feet. Default 0.
            is_map_centered_by_origin: Determines if the map should be centered at the pre-offset origin. Default False.
            label: A custom label string for the file.
            zoom: The zoom level of the map in Google Maps API.
            should_confirm: If the method should prompt the user to confirm map generation, to ensure non-accidental usage of the Google Maps API. Default True.
        Returns:
            A file as a jpg with the radius using specified settings. Name is formatted as label + radius + latOffset + longOffset.
        """
        # Confirm that the function wasn't called by accident and avoid billing hike
        if should_confirm:
            confirm = input("Do you want to generate a map with the title \"" + label + "_" + str("TBD" if radius == None else radius) + "_" + str(latOffset) + "_" + str(longOffset) + ".jpg" + "\"? (y/N) ").lower()
            if not confirm == "y" and not confirm == "yes": return

        # proceed with the rest of the function
        origin = self.origin if input_origin == None else input_origin
        offsetOrigin = self.returnOffset(origin, latFeet=latOffset, longFeet=longOffset) if (latOffset != 0 or longOffset != 0) else origin
        circular_points = self.get_circular_points(radius=radius, origin=offsetOrigin, num=num, phi_offset=(360/num)*phi_offset)
        print(len(circular_points))
        path_points = self.get_elevation_dict(circular_points, offsetOrigin, radius, division)
        print(path_points)
        path = maps.StaticMapPath(points=self.get_path_points_as_array(path_points), weight=3)
        waypoint = maps.StaticMapMarker(maps.convert.latlng(offsetOrigin), 5)
        mappa = maps.static_map(self.key_client, (2000, 2000), maps.convert.latlng(origin if is_map_centered_by_origin else offsetOrigin), maptype="satellite", path=path, zoom=zoom, markers=waypoint)
        print(type(mappa))
        avg_radius = [rad["radius"] for rad in circular_points.values()]
        avg_radius = sum(avg_radius) / len(avg_radius)
        f = open(label + "_" + str(avg_radius) + "_" + str(latOffset) + "_" + str(longOffset) + "_" + str(phi_offset) + "_" + ".jpg", 'wb')
        for chunk in mappa:
            if chunk:
                f.write(chunk)
        f.close()

# Testing/Usage
fragmentation = HFD()

r, m = fragmentation.get_radius_and_m(
    ethane_tank_mass=22.2,
    ethane_mass=4.99,
    nitrous_tank_mass=7.0,
    nitrous_mass=2.72,
    yield_factor=1.0,
    fuel_tank_diameter=0.18,
    fuel_tank_length=0.84,
    fuel_tank_density=2.710,
    oxidizer_tank_diameter=0.18,
    oxidizer_tank_length=0.41,
    oxidizer_tank_density=2.710,
    only_tanks=False)


r_small, m_small = fragmentation.get_radius_and_m(
    ethane_tank_mass=7.0,#22.2,
    ethane_mass=2.05,#4.99,
    nitrous_tank_mass=7.0,
    nitrous_mass=2.72,
    yield_factor=1.0,
    fuel_tank_diameter=0.18,
    fuel_tank_length=0.41,#0.84,
    fuel_tank_density=2.710,
    oxidizer_tank_diameter=0.18,
    oxidizer_tank_length=0.41,
    oxidizer_tank_density=2.710,
    only_tanks=False)


rocket_r, rocket_m = fragmentation.get_radius_and_m(
    ethane_tank_mass=2.329806077941931,
    ethane_mass=6.7,
    nitrous_tank_mass=4.97,
    nitrous_mass=1.3064,
    yield_factor=1.0,
    fuel_tank_diameter=0.15,
    fuel_tank_length=0.3,
    fuel_tank_density=2.710,
    oxidizer_tank_diameter=0.15,
    oxidizer_tank_length=0.64,
    oxidizer_tank_density=2.710,
    only_tanks=False)

gse_r, gse_m =  fragmentation.get_radius_and_m(
    ethane_tank_mass=54.0,
    ethane_mass=14.51,
    nitrous_tank_mass=54.0,
    nitrous_mass=27.2,
    yield_factor=1.0,
    fuel_tank_diameter=0.23,
    fuel_tank_length=1.30,
    fuel_tank_density=7.75,
    oxidizer_tank_diameter=0.23,
    oxidizer_tank_length=1.30,
    oxidizer_tank_density=7.75,
    only_tanks=False)

blast_mapper = Blast()
blast_mapper.map_points(radius=None, longOffset=0.0, latOffset=0.0, num=30, label="stony_creek_hdf_fragement_calculations_ehs_draft_phi", zoom=19, phi_offset=0)