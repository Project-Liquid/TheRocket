import os
import sys
import csv
import threading
import math
from collections import deque
from datetime import datetime

import serial
import serial.tools.list_ports


from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QComboBox, QSpinBox,
    QGroupBox, QSizePolicy, QFileDialog, QMessageBox, QFrame, QTextEdit, QLineEdit
)
from PyQt5.QtCore import Qt, QTimer, pyqtSignal, QObject, QThread
from PyQt5.QtGui import QFont, QColor, QPalette

import pyqtgraph as pg

import pickle
from pathlib import Path


# ─────────────────────────────────────────
#  Serial Worker (runs in background thread)
# ─────────────────────────────────────────
class SerialWorker(QObject):
    data_received = pyqtSignal(dict)
    raw_received = pyqtSignal(str)  # Raw serial line for monitor
    connection_lost = pyqtSignal()

    def __init__(self, port, baud=57600): # New baud rate
        super().__init__()
        self.port = port
        self.baud = baud
        self._running = False
        self.ser = None

    def start(self):
        self._running = True
        try:
            self.ser = serial.Serial(self.port, self.baud, timeout=1)
        except serial.SerialException as e:
            print(f"Serial error: {e}")
            self.connection_lost.emit()
            return

        num_raw = 0
        while self._running:
            try:
                raw = self.ser.readline().decode('ascii', errors='replace').rstrip('\r\n')
                print(raw)
                if raw.startswith("DATA|"):
                    parsed = self._parse(raw)
                    if parsed:
                        self.data_received.emit(parsed)
                    self.raw_received.emit(raw)  # Emit raw line for monitor
                elif raw != "":  # Emit non-empty lines that don't start with DATA| as well
                    self.raw_received.emit(raw)
                    


            except (serial.SerialException, OSError):
                self.connection_lost.emit()
                break

    def send(self, cmd: str):
        if self.ser and self.ser.is_open:
            self.ser.write((cmd + '\n').encode())

    def stop(self):
        self._running = False
        if self.ser and self.ser.is_open:
            self.ser.close()

    @staticmethod
    def _parse(line: str) -> dict:
        """Parse DATA|elapsed_s|KEY:VAL|KEY:VAL|... into a dict.

        The Arduino sends elapsed seconds as a float (e.g. 43.52).
        PT values may be empty strings if the firmware omits Serial.print()
        before the delimiter — those keys are skipped gracefully.
        """
        try:
            parts = line.split('|')
            # parts[0] = "DATA", parts[1] = elapsed seconds (float)
            result = {'millis': float(parts[1]) * 1000}  # keep as ms for compat
            for part in parts[2:]:
                if ':' in part:
                    k, v = part.split(':', 1)
                    v = v.strip()
                    if v == '':        # firmware didn't emit a value — skip key
                        continue
                    try:
                        result[k] = float(v)
                    except ValueError: # non-numeric token — skip
                        continue
            return result
        except Exception:
            return {}


# ─────────────────────────────────────────
#  Valve Button Widget
# ─────────────────────────────────────────
class ValveButton(QPushButton):
    def __init__(self, label, cmd_on, cmd_off, parent=None):
        super().__init__(label, parent)
        self.cmd_on = cmd_on
        self.cmd_off = cmd_off
        self._open = False
        self.setFixedHeight(48)
        self.setFont(QFont("Courier New", 10, QFont.Bold))
        self._refresh()

    def toggle(self, send_fn):
        self._open = not self._open
        send_fn(self.cmd_on if self._open else self.cmd_off)
        self._refresh()

    def set_state(self, is_open: bool):
        self._open = is_open
        self._refresh()

    def _refresh(self):
        if self._open:
            self.setStyleSheet(
                "background:#00ff88; color:#000; border:2px solid #00cc66;"
                "border-radius:4px; font-weight:bold;"
            )
            self.setText(self.text().split('●')[0].strip() + "  ● OPEN")
        else:
            self.setStyleSheet(
                "background:#1a1a2e; color:#ff4466; border:2px solid #ff4466;"
                "border-radius:4px; font-weight:bold;"
            )
            self.setText(self.text().split('●')[0].strip().replace('  ', '') + "  ● CLOSED")


# ─────────────────────────────────────────
#  Sensor Readout Widget
# ─────────────────────────────────────────
class SensorLabel(QFrame):
    def __init__(self, title, unit, warning_hi=None, parent=None):
        super().__init__(parent)
        self.unit = unit
        self.warning_hi = warning_hi
        self.setFrameShape(QFrame.StyledPanel)
        self.setStyleSheet("background:#0d1117; border:1px solid #30363d; border-radius:6px;")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 6, 8, 6)

        self.title_lbl = QLabel(title)
        self.title_lbl.setFont(QFont("Courier New", 8))
        self.title_lbl.setStyleSheet("color:#8b949e; border:none;")

        self.value_lbl = QLabel("---")
        self.value_lbl.setFont(QFont("Courier New", 18, QFont.Bold))
        self.value_lbl.setStyleSheet("color:#58a6ff; border:none;")
        self.value_lbl.setAlignment(Qt.AlignRight)

        self.unit_lbl = QLabel(unit)
        self.unit_lbl.setFont(QFont("Courier New", 8))
        self.unit_lbl.setStyleSheet("color:#8b949e; border:none;")
        self.unit_lbl.setAlignment(Qt.AlignRight)

        layout.addWidget(self.title_lbl)
        layout.addWidget(self.value_lbl)
        layout.addWidget(self.unit_lbl)

    def update_value(self, val: float):
        self.value_lbl.setText(f"{val:.1f}")
        if self.warning_hi and val > self.warning_hi:
            self.value_lbl.setStyleSheet("color:#ff4466; border:none;")
        else:
            self.value_lbl.setStyleSheet("color:#58a6ff; border:none;")


# ─────────────────────────────────────────
#  Main Window
# ─────────────────────────────────────────
class GroundStation(QMainWindow):
    HISTORY_LEN = 300  # samples kept for chart (~150s at 500ms)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("GSE Ground Station")
        self.setMinimumSize(1280, 800)
        self._apply_dark_theme()

        self.worker = None
        self.serial_thread = None
        self.send_fn = lambda cmd: None  # no-op until connected

        # Data history for charts
        self.t_hist    = deque(maxlen=self.HISTORY_LEN)
        self.et_up_hist = deque(maxlen=self.HISTORY_LEN)
        self.et_dn_hist = deque(maxlen=self.HISTORY_LEN)
        self.nit_up_hist = deque(maxlen=self.HISTORY_LEN)
        self.nit_dn_hist = deque(maxlen=self.HISTORY_LEN)
        self.lc_et_hist = deque(maxlen=self.HISTORY_LEN)
        self.lc_nit_hist = deque(maxlen=self.HISTORY_LEN)

        # CSV logging
        self.log_rows = []
        self.logging_active = False

        # Serial monitor
        self.serial_monitor_lines = deque(maxlen=500)

        # MBV positions
        self.mbv_e_pos = 0.0
        self.mbv_n_pos = 0.0

        self._build_ui()

        # Folder directory to save data
        self.parent_folder = Path(Path.cwd().as_posix() + "/stream_data")
        self.log_start_time: datetime | None = None
        self.log_folder_path: Path | None = None


    # ── UI Construction ──────────────────────────────────────────
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setSpacing(8)
        root.setContentsMargins(12, 12, 12, 12)

        # Top bar: connection
        root.addLayout(self._build_connection_bar())

        # Main content
        content = QHBoxLayout()
        content.setSpacing(8)

        left = QVBoxLayout()
        left.addWidget(self._build_pt_group())
        left.addWidget(self._build_lc_group())
        left.addWidget(self._build_serial_monitor_group())
        left.addWidget(self._build_export_group())
        left.addStretch()

        right = QVBoxLayout()
        right.addWidget(self._build_valve_group())
        right.addWidget(self._build_mbv_group())
        right.addWidget(self._build_cold_flow_group())
        right.addWidget(self._build_static_fire_group())
        right.addStretch()

        content.addLayout(left, 3)
        content.addWidget(self._build_chart_group(), 5)
        content.addLayout(right, 2)

        root.addLayout(content)

    def _build_connection_bar(self):
        bar = QHBoxLayout()

        lbl = QLabel("PORT")
        lbl.setFont(QFont("Courier New", 9))
        lbl.setStyleSheet("color:#8b949e;")

        self.port_combo = QComboBox()
        self.port_combo.setFixedWidth(160)
        self._refresh_ports()

        refresh_btn = QPushButton("↺ Refresh")
        refresh_btn.setFixedWidth(90)
        refresh_btn.clicked.connect(self._refresh_ports)

        self.connect_btn = QPushButton("Connect")
        self.connect_btn.setFixedWidth(100)
        self.connect_btn.clicked.connect(self._toggle_connection)
        self.connect_btn.setStyleSheet("background:#238636; color:#fff; border-radius:4px;")

        self.status_lbl = QLabel("● Disconnected")
        self.status_lbl.setFont(QFont("Courier New", 9, QFont.Bold))
        self.status_lbl.setStyleSheet("color:#ff4466;")

        bar.addWidget(lbl)
        bar.addWidget(self.port_combo)
        bar.addWidget(refresh_btn)
        bar.addWidget(self.connect_btn)
        bar.addWidget(self.status_lbl)
        bar.addStretch()
        return bar

    def _build_pt_group(self):
        grp = QGroupBox("PRESSURE TRANSDUCERS")
        grp.setFont(QFont("Courier New", 9, QFont.Bold))
        grid = QGridLayout(grp)
        grid.setSpacing(6)

        # TODO: The redlines in GUI are not the same as redline in GUI. Must change
        # once redlines are determined by fluids.
        self.pt_et_up  = SensorLabel("Ethane Upstream",   "psi", warning_hi=900)
        self.pt_et_dn  = SensorLabel("Ethane Downstream", "psi", warning_hi=900)
        self.pt_nit_up = SensorLabel("Nitrous Upstream",  "psi", warning_hi=1400)
        self.pt_nit_dn = SensorLabel("Nitrous Downstream","psi", warning_hi=1400)

        grid.addWidget(self.pt_et_up,  0, 0)
        grid.addWidget(self.pt_et_dn,  0, 1)
        grid.addWidget(self.pt_nit_up, 1, 0)
        grid.addWidget(self.pt_nit_dn, 1, 1)
        return grp

    def _build_lc_group(self):
        grp = QGroupBox("LOAD CELLS")
        grp.setFont(QFont("Courier New", 9, QFont.Bold))
        grid = QGridLayout(grp)
        grid.setSpacing(6)

        self.lc_e1 = SensorLabel("Ethane LC1",  "lbs")
        self.lc_e2 = SensorLabel("Ethane LC2",  "lbs")
        self.lc_e3 = SensorLabel("Ethane LC3",  "lbs")
        self.lc_et = SensorLabel("Ethane TOTAL","lbs")
        self.lc_n1 = SensorLabel("Nitrous LC1", "lbs")
        self.lc_n2 = SensorLabel("Nitrous LC2", "lbs")
        self.lc_n3 = SensorLabel("Nitrous LC3", "lbs")
        self.lc_nt = SensorLabel("Nitrous TOTAL","lbs")
        self.lc_t = SensorLabel("Thrust LC", "lbs")

        for i, w in enumerate([self.lc_e1, self.lc_e2, self.lc_e3, self.lc_et]):
            grid.addWidget(w, 0, i)
        for i, w in enumerate([self.lc_n1, self.lc_n2, self.lc_n3, self.lc_nt]):
            grid.addWidget(w, 1, i)
        grid.addWidget(self.lc_t, 2, 1, 1, 2)
        return grp

    def _build_valve_group(self):
        grp = QGroupBox("VALVES")
        grp.setFont(QFont("Courier New", 9, QFont.Bold))
        layout = QVBoxLayout(grp)
        layout.setSpacing(6)

        self.btn_erv = ValveButton("Ethane Run Valve",   "ETHANE_RUN_ON",   "ETHANE_RUN_OFF")
        self.btn_ev  = ValveButton("Ethane Vent",        "ETHANE_VENT_ON",  "ETHANE_VENT_OFF")
        self.btn_nrv = ValveButton("Nitrous Run Valve",  "NITROUS_RUN_ON",  "NITROUS_RUN_OFF")
        self.btn_nv  = ValveButton("Nitrous Vent",       "NITROUS_VENT_ON", "NITROUS_VENT_OFF")

        for btn in [self.btn_erv, self.btn_ev, self.btn_nrv, self.btn_nv]:
            btn.clicked.connect(lambda checked, b=btn: b.toggle(self.send_fn))
            layout.addWidget(btn)

        # Emergency stop
        estop = QPushButton("⚠  EMERGENCY STOP")
        estop.setFixedHeight(52)
        estop.setFont(QFont("Courier New", 11, QFont.Bold))
        estop.setStyleSheet(
            "background:#ff4466; color:#fff; border-radius:4px; border:2px solid #ff0033;"
        )
        estop.clicked.connect(self._emergency_stop)
        layout.addWidget(estop)

        return grp

    def _build_mbv_group(self):
        grp = QGroupBox("MOTORIZED BALL VALVES")
        grp.setFont(QFont("Courier New", 9, QFont.Bold))
        layout = QVBoxLayout(grp)
        layout.setSpacing(6)

        # Ethane MBV
        ethane_layout = QVBoxLayout()
        ethane_lbl = QLabel("Ethane MBV")
        ethane_lbl.setFont(QFont("Courier New", 9, QFont.Bold))
        ethane_lbl.setStyleSheet("color:#58a6ff;")
        self.mbv_e_display = SensorLabel("Position", "°")
        ethane_layout.addWidget(ethane_lbl)
        ethane_layout.addWidget(self.mbv_e_display)
        
        ethane_btn_layout = QHBoxLayout()
        btn_e_10 = QPushButton("+10°")
        btn_e_10.setFont(QFont("Courier New", 9))
        btn_e_10.clicked.connect(lambda: self.send_fn("ETHANE_MBV_10"))
        btn_e_90 = QPushButton("+90°")
        btn_e_90.setFont(QFont("Courier New", 9))
        btn_e_90.clicked.connect(lambda: self.send_fn("ETHANE_MBV_90"))
        ethane_btn_layout.addWidget(btn_e_10)
        ethane_btn_layout.addWidget(btn_e_90)
        ethane_layout.addLayout(ethane_btn_layout)

        # Nitrous MBV
        nitrous_layout = QVBoxLayout()
        nitrous_lbl = QLabel("Nitrous MBV")
        nitrous_lbl.setFont(QFont("Courier New", 9, QFont.Bold))
        nitrous_lbl.setStyleSheet("color:#ff7b72;")
        self.mbv_n_display = SensorLabel("Position", "°")
        nitrous_layout.addWidget(nitrous_lbl)
        nitrous_layout.addWidget(self.mbv_n_display)
        
        nitrous_btn_layout = QHBoxLayout()
        btn_n_10 = QPushButton("+10°")
        btn_n_10.setFont(QFont("Courier New", 9))
        btn_n_10.clicked.connect(lambda: self.send_fn("NITROUS_MBV_10"))
        btn_n_90 = QPushButton("+90°")
        btn_n_90.setFont(QFont("Courier New", 9))
        btn_n_90.clicked.connect(lambda: self.send_fn("NITROUS_90"))
        nitrous_btn_layout.addWidget(btn_n_10)
        nitrous_btn_layout.addWidget(btn_n_90)
        nitrous_layout.addLayout(nitrous_btn_layout)

        layout.addLayout(ethane_layout)
        layout.addLayout(nitrous_layout)
        return grp

    def _build_cold_flow_group(self):
        grp = QGroupBox("COLD FLOW AUTOMATED")
        grp.setFont(QFont("Courier New", 9, QFont.Bold))
        layout = QVBoxLayout(grp)
        layout.setSpacing(6)

        side_layout = QHBoxLayout()
        side_lbl = QLabel("Side:")
        side_lbl.setFont(QFont("Courier New", 9))
        self.cf_side = QComboBox()
        self.cf_side.addItems(["ETHANE", "NITROUS"])
        side_layout.addWidget(side_lbl)
        side_layout.addWidget(self.cf_side)

        time_layout = QHBoxLayout()
        time_lbl = QLabel("Run time (s):")
        time_lbl.setFont(QFont("Courier New", 9))
        self.cf_time = QSpinBox()
        self.cf_time.setRange(1, 600)
        self.cf_time.setValue(10)
        time_layout.addWidget(time_lbl)
        time_layout.addWidget(self.cf_time)

        self.cf_btn = QPushButton("▶  START COLD FLOW")
        self.cf_btn.setFixedHeight(42)
        self.cf_btn.setFont(QFont("Courier New", 10, QFont.Bold))
        self.cf_btn.setStyleSheet(
            "background:#1f6feb; color:#fff; border-radius:4px; border:2px solid #388bfd;"
        )
        self.cf_btn.clicked.connect(self._start_cold_flow)

        self.cf_status = QLabel("Ready")
        self.cf_status.setFont(QFont("Courier New", 8))
        self.cf_status.setStyleSheet("color:#8b949e;")

        layout.addLayout(side_layout)
        layout.addLayout(time_layout)
        layout.addWidget(self.cf_btn)
        layout.addWidget(self.cf_status)
        return grp

    def _build_static_fire_group(self):
        grp = QGroupBox("STATIC FIRE")
        grp.setFont(QFont("Courier New", 9, QFont.Bold))
        layout = QVBoxLayout(grp)
        layout.setSpacing(6)

        time_layout = QHBoxLayout()
        time_lbl = QLabel("Run time (s):")
        time_lbl.setFont(QFont("Courier New", 9))
        self.sf_time = QSpinBox()
        self.sf_time.setRange(1, 600)
        self.sf_time.setValue(10)
        time_layout.addWidget(time_lbl)
        time_layout.addWidget(self.sf_time)

        self.sf_btn = QPushButton("🔥  START STATIC FIRE")
        self.sf_btn.setFixedHeight(42)
        self.sf_btn.setFont(QFont("Courier New", 10, QFont.Bold))
        self.sf_btn.setStyleSheet(
            "background:#1f6feb; color:#fff; border-radius:4px; border:2px solid #388bfd;"
        )
        self.sf_btn.clicked.connect(self._start_static_fire)

        self.sf_status = QLabel("Ready")
        self.sf_status.setFont(QFont("Courier New", 8))
        self.sf_status.setStyleSheet("color:#8b949e;")

        layout.addLayout(time_layout)
        layout.addWidget(self.sf_btn)
        layout.addWidget(self.sf_status)
        return grp

    def _build_export_group(self):
        grp = QGroupBox("DATA LOGGING")
        grp.setFont(QFont("Courier New", 9, QFont.Bold))
        layout = QVBoxLayout(grp)

        self.log_btn = QPushButton("⏺  Start Logging")
        self.log_btn.setFixedHeight(36)
        self.log_btn.setFont(QFont("Courier New", 9, QFont.Bold))
        self.log_btn.setStyleSheet("background:#388bfd; color:#fff; border-radius:4px;")
        self.log_btn.clicked.connect(self._toggle_logging)

        export_btn = QPushButton("💾  Export CSV")
        export_btn.setFixedHeight(36)
        export_btn.setFont(QFont("Courier New", 9, QFont.Bold))
        export_btn.setStyleSheet("background:#238636; color:#fff; border-radius:4px;")
        export_btn.clicked.connect(self._export_csv)

        load_backup_btn = QPushButton("📂  Load Backup")
        load_backup_btn.setFixedHeight(36)
        load_backup_btn.setFont(QFont("Courier New", 9, QFont.Bold))
        load_backup_btn.setStyleSheet("background:#8b949e; color:#fff; border-radius:4px;")
        load_backup_btn.clicked.connect(self._load_backup)

        self.log_count_lbl = QLabel("0 rows logged")
        self.log_count_lbl.setFont(QFont("Courier New", 8))
        self.log_count_lbl.setStyleSheet("color:#8b949e;")

        layout.addWidget(self.log_btn)
        layout.addWidget(export_btn)
        layout.addWidget(load_backup_btn)
        layout.addWidget(self.log_count_lbl)
        return grp

    def _build_serial_monitor_group(self):
        grp = QGroupBox("SERIAL MONITOR")
        grp.setFont(QFont("Courier New", 9, QFont.Bold))
        layout = QVBoxLayout(grp)

        self.serial_monitor = QTextEdit()
        self.serial_monitor.setReadOnly(True)
        self.serial_monitor.setFont(QFont("Courier New", 8))
        self.serial_monitor.setStyleSheet(
            "background:#0d1117; color:#58a6ff; border:1px solid #30363d; border-radius:4px;"
        )
        self.serial_monitor.setMinimumHeight(150)

        # Clear button
        clear_btn = QPushButton("Clear")
        clear_btn.setFixedHeight(28)
        clear_btn.setFont(QFont("Courier New", 8))
        clear_btn.clicked.connect(self.serial_monitor.clear)

        # Input section for sending commands
        input_layout = QHBoxLayout()
        input_layout.setSpacing(4)
        
        self.serial_input = QLineEdit()
        self.serial_input.setFont(QFont("Courier New", 8))
        self.serial_input.setStyleSheet(
            "background:#161b22; color:#c9d1d9; border:1px solid #30363d; border-radius:4px; padding:4px;"
        )
        self.serial_input.setPlaceholderText("Enter command...")
        self.serial_input.returnPressed.connect(self._send_serial_command)
        
        send_btn = QPushButton("Send")
        send_btn.setFixedHeight(28)
        send_btn.setFixedWidth(60)
        send_btn.setFont(QFont("Courier New", 8))
        send_btn.setStyleSheet("background:#238636; color:#fff; border-radius:4px;")
        send_btn.clicked.connect(self._send_serial_command)
        
        input_layout.addWidget(self.serial_input)
        input_layout.addWidget(send_btn)

        layout.addWidget(self.serial_monitor)
        layout.addWidget(clear_btn)
        layout.addLayout(input_layout)
        return grp

    def _build_chart_group(self):
        grp = QGroupBox("LIVE DATA")
        grp.setFont(QFont("Courier New", 9, QFont.Bold))
        layout = QVBoxLayout(grp)

        pg.setConfigOption('background', '#0d1117')
        pg.setConfigOption('foreground', '#8b949e')

        self.pt_chart = pg.PlotWidget(title="Pressure (psi)")
        self.pt_chart.addLegend()
        self.pt_chart.showGrid(x=True, y=True, alpha=0.2)
        self.curve_et_up  = self.pt_chart.plot(pen=pg.mkPen('#58a6ff', width=2), name="Eth Up")
        self.curve_et_dn  = self.pt_chart.plot(pen=pg.mkPen('#79c0ff', width=1), name="Eth Dn")
        self.curve_nit_up = self.pt_chart.plot(pen=pg.mkPen('#ff7b72', width=2), name="Nit Up")
        self.curve_nit_dn = self.pt_chart.plot(pen=pg.mkPen('#ffa198', width=1), name="Nit Dn")

        self.lc_chart = pg.PlotWidget(title="Propellant Mass (lbs)")
        self.lc_chart.addLegend()
        self.lc_chart.showGrid(x=True, y=True, alpha=0.2)
        self.curve_lc_et  = self.lc_chart.plot(pen=pg.mkPen('#3fb950', width=2), name="Ethane")
        self.curve_lc_nit = self.lc_chart.plot(pen=pg.mkPen('#d29922', width=2), name="Nitrous")

        clear_charts_btn = QPushButton("Clear Charts")
        clear_charts_btn.setFixedHeight(32)
        clear_charts_btn.setFont(QFont("Courier New", 8))
        clear_charts_btn.setStyleSheet("background:#da3633; color:#fff; border-radius:4px;")
        clear_charts_btn.clicked.connect(self._clear_chart_display)

        layout.addWidget(self.pt_chart)
        layout.addWidget(self.lc_chart)
        layout.addWidget(clear_charts_btn)
        return grp

    # ── Connection ────────────────────────────────────────────────
    def _refresh_ports(self):
        self.port_combo.clear()
        ports = [p.device for p in serial.tools.list_ports.comports()]
        self.port_combo.addItems(ports if ports else ["No ports found"])

    def _toggle_connection(self):
        if self.serial_thread and self.serial_thread.isRunning():
            self._disconnect()
        else:
            self._connect()

    def _connect(self):
        port = self.port_combo.currentText()
        if not port or port == "No ports found":
            QMessageBox.warning(self, "No Port", "Select a valid serial port.")
            return

        self.worker = SerialWorker(port)
        self.serial_thread = QThread()
        self.worker.moveToThread(self.serial_thread)
        self.serial_thread.started.connect(self.worker.start)
        self.worker.data_received.connect(self._on_data)
        self.worker.raw_received.connect(self._on_raw_data)
        self.worker.connection_lost.connect(self._disconnect)
        self.send_fn = self.worker.send
        self.serial_thread.start()

        self.connect_btn.setText("Disconnect")
        self.connect_btn.setStyleSheet("background:#da3633; color:#fff; border-radius:4px;")
        self.status_lbl.setText("● Connected")
        self.status_lbl.setStyleSheet("color:#3fb950;")


    def _backup_connect(self):
        port = self.port_combo.currentText()
        if not port or port == "No ports found":
            QMessageBox.warning(self, "No Port", "Select a valid serial port.")
            return

        self.worker = SerialWorker(port)
        self.serial_thread = QThread()
        self.worker.moveToThread(self.serial_thread)
        self.serial_thread.started.connect(self.worker.start)
        self.worker.data_received.connect(self._on_data)
        self.worker.raw_received.connect(self._on_raw_data)
        self.worker.connection_lost.connect(self._disconnect)
        self.send_fn = self.worker.send
        self.serial_thread.start()

        self.connect_btn.setText("Close Historical View")
        self.connect_btn.setStyleSheet("background:#da3633; color:#fff; border-radius:4px;")
        self.status_lbl.setText("● Backup Loaded")
        self.status_lbl.setStyleSheet("color:#3fb950;")


    # TODO: If CSV is recording, make sure to finalize and save the file on disconnect
    def _disconnect(self):
        if self.worker:
            self.worker.stop()
        if self.serial_thread:
            self.serial_thread.quit()
            self.serial_thread.wait()
        self.send_fn = lambda cmd: None
        self.connect_btn.setText("Connect")
        self.connect_btn.setStyleSheet("background:#238636; color:#fff; border-radius:4px;")
        self.status_lbl.setText("● Disconnected")
        self.status_lbl.setStyleSheet("color:#ff4466;")
        self.serial_monitor.clear()

    # ── Data Handler ──────────────────────────────────────────────
    def _on_data(self, state: dict, packet_size=100):
        t = state.get('millis', 0) / 1000.0

        # PT readouts
        et_up  = state.get('PT_EU', float('nan'))
        et_dn  = state.get('PT_ED', float('nan'))
        nit_up = state.get('PT_NU', float('nan'))
        nit_dn = state.get('PT_ND', float('nan'))
        if 'PT_EU' in state:
            self.pt_et_up.update_value(et_up)
        if 'PT_ED' in state:
            self.pt_et_dn.update_value(et_dn)
        if 'PT_NU' in state:
            self.pt_nit_up.update_value(nit_up)
        if 'PT_ND' in state:
            self.pt_nit_dn.update_value(nit_dn)

        # LC readouts
        lc1 = state.get('LC_E1', 0.0) if 'LC_E1' in state else 0.0
        lc2 = state.get('LC_E2', 0.0) if 'LC_E2' in state else 0.0
        lc3 = state.get('LC_E3', 0.0) if 'LC_E3' in state else 0.0
        et_total = lc1 + lc2 + lc3 if all(k in state for k in ['LC_E1', 'LC_E2', 'LC_E3']) else float('nan')
        n1  = state.get('LC_N1', 0.0) if 'LC_N1' in state else 0.0
        n2  = state.get('LC_N2', 0.0) if 'LC_N2' in state else 0.0
        n3  = state.get('LC_N3', 0.0) if 'LC_N3' in state else 0.0
        nit_total = n1 + n2 + n3 if all(k in state for k in ['LC_N1', 'LC_N2', 'LC_N3']) else float('nan')
        
        if 'LC_E1' in state:
            self.lc_e1.update_value(lc1)
        if 'LC_E2' in state:
            self.lc_e2.update_value(lc2)
        if 'LC_E3' in state:
            self.lc_e3.update_value(lc3)
        if all(k in state for k in ['LC_E1', 'LC_E2', 'LC_E3']):
            self.lc_et.update_value(et_total)
        if 'LC_N1' in state:
            self.lc_n1.update_value(n1)
        if 'LC_N2' in state:
            self.lc_n2.update_value(n2)
        if 'LC_N3' in state:
            self.lc_n3.update_value(n3)
        if all(k in state for k in ['LC_N1', 'LC_N2', 'LC_N3']):
            self.lc_nt.update_value(nit_total)

        # Thrust load cell
        lc_t = state.get('LC_T', float('nan'))
        if 'LC_T' in state:
            self.lc_t.update_value(lc_t)

        # MBV positions
        mbv_e = state.get('MBV_E', float('nan'))
        mbv_n = state.get('MBV_N', float('nan'))
        self.mbv_e_pos = mbv_e
        self.mbv_n_pos = mbv_n
        if 'MBV_E' in state:
            self.mbv_e_display.update_value(mbv_e)
        if 'MBV_N' in state:
            self.mbv_n_display.update_value(mbv_n)

        # Valve states from firmware flags (re-enable if firmware emits these)
        if 'ERV' in state:
            self.btn_erv.set_state(bool(state.get('ERV', 0)))
        if 'EV' in state:
            self.btn_ev.set_state(bool(state.get('EV', 0)))
        if 'NRV' in state:
            self.btn_nrv.set_state(bool(state.get('NRV', 0)))
        if 'NV' in state:
            self.btn_nv.set_state(bool(state.get('NV', 0)))

        # Charts
        self.t_hist.append(t)
        self.et_up_hist.append(et_up)
        self.et_dn_hist.append(et_dn)
        self.nit_up_hist.append(nit_up)
        self.nit_dn_hist.append(nit_dn)
        self.lc_et_hist.append(et_total)
        self.lc_nit_hist.append(nit_total)

        # Only plot valid (non-NaN) values
        tl = list(self.t_hist)
        
        # Filter NaN values for each curve
        def plot_valid(t_list, val_list):
            valid_t = [t for t, v in zip(t_list, val_list) if not math.isnan(v)]
            valid_v = [v for v in val_list if not math.isnan(v)]
            return valid_t, valid_v
        
        t_et_up, v_et_up = plot_valid(tl, list(self.et_up_hist))
        t_et_dn, v_et_dn = plot_valid(tl, list(self.et_dn_hist))
        t_nit_up, v_nit_up = plot_valid(tl, list(self.nit_up_hist))
        t_nit_dn, v_nit_dn = plot_valid(tl, list(self.nit_dn_hist))
        t_lc_et, v_lc_et = plot_valid(tl, list(self.lc_et_hist))
        t_lc_nit, v_lc_nit = plot_valid(tl, list(self.lc_nit_hist))
        
        self.curve_et_up.setData(t_et_up, v_et_up)
        self.curve_et_dn.setData(t_et_dn, v_et_dn)
        self.curve_nit_up.setData(t_nit_up, v_nit_up)
        self.curve_nit_dn.setData(t_nit_dn, v_nit_dn)
        self.curve_lc_et.setData(t_lc_et, v_lc_et)
        self.curve_lc_nit.setData(t_lc_nit, v_lc_nit)

        # CSV logging
        if self.logging_active:
            row = {'time_s': t}
            if 'PT_EU' in state:
                row['ET_UP'] = state['PT_EU']
            if 'PT_ED' in state:
                row['ET_DN'] = state['PT_ED']
            if 'PT_NU' in state:
                row['NIT_UP'] = state['PT_NU']
            if 'PT_ND' in state:
                row['NIT_DN'] = state['PT_ND']
            if 'LC_E1' in state:
                row['LC1'] = state['LC_E1']
            if 'LC_E2' in state:
                row['LC2'] = state['LC_E2']
            if 'LC_E3' in state:
                row['LC3'] = state['LC_E3']
            if all(k in state for k in ['LC_E1', 'LC_E2', 'LC_E3']):
                row['ET_TOTAL'] = et_total
            if 'LC_N1' in state:
                row['NLC1'] = state['LC_N1']
            if 'LC_N2' in state:
                row['NLC2'] = state['LC_N2']
            if 'LC_N3' in state:
                row['NLC3'] = state['LC_N3']
            if all(k in state for k in ['LC_N1', 'LC_N2', 'LC_N3']):
                row['NIT_TOTAL'] = nit_total
            if 'LC_T' in state:
                row['LC_T'] = state['LC_T']
            if 'ERV' in state:
                row['ERV'] = state['ERV']
            if 'EV' in state:
                row['EV'] = state['EV']
            if 'NRV' in state:
                row['NRV'] = state['NRV']
            if 'NV' in state:
                row['NV'] = state['NV']
            self.log_rows.append(row)

            # Backup Data Storage
            num_rows = len(self.log_rows)
            self.log_count_lbl.setText(f"{num_rows} rows logged")
            # PACKET_SIZE = 100 # how many rows to save per pickle file

            if self.logging_active and num_rows % packet_size == 0 and self.log_folder_path is not None:
                path = self.log_folder_path / f"data-{num_rows}.pkl"

                # Create folder if it doesn't exist
                path.parent.mkdir(parents=True, exist_ok=True)
                with open(path, "wb") as f:
                    pickle.dump(self.log_rows[-packet_size:], f)


    def _on_raw_data(self, raw: str):
        """Append raw serial line to monitor."""
        self.serial_monitor.append(raw)

        if self.logging_active and not raw.startswith("DATA|") and self.log_folder_path is not None and self.log_rows is not None:  # skip echo lines
            path = self.log_folder_path / f"raw-{len(self.log_rows)}.pkl"

            # Create folder if it doesn't exist
            path.parent.mkdir(parents=True, exist_ok=True)
            with open(path, "wb") as f:
                pickle.dump(raw, f)
            
        # Auto-scroll to bottom
        self.serial_monitor.verticalScrollBar().setValue(
            self.serial_monitor.verticalScrollBar().maximum()
        )

    # ── Controls ──────────────────────────────────────────────────
    def _emergency_stop(self):
        for cmd in ["ETHANE_RUN_OFF", "NITROUS_RUN_OFF",
                    "ETHANE_VENT_ON", "NITROUS_VENT_ON"]:
            self.send_fn(cmd)

    def _send_serial_command(self):
        """Send command from serial monitor input field."""
        cmd = self.serial_input.text().strip()
        if cmd:
            self.send_fn(cmd)
            self.serial_monitor.append(f"→ {cmd}")
            self.serial_input.clear()

    def _clear_chart_display(self):
        """Clear chart display without erasing historical data."""
        self.curve_et_up.setData([], [])
        self.curve_et_dn.setData([], [])
        self.curve_nit_up.setData([], [])
        self.curve_nit_dn.setData([], [])
        self.curve_lc_et.setData([], [])
        self.curve_lc_nit.setData([], [])

    def _start_cold_flow(self):
        side     = self.cf_side.currentText()
        run_time = self.cf_time.value()
        reply = QMessageBox.question(
            self, "Confirm Cold Flow",
            f"Run cold flow on {side} side for {run_time}s?\n\nThis will open valves.",
            QMessageBox.Yes | QMessageBox.No
        )
        if reply == QMessageBox.Yes:
            run_time_ms = run_time * 1000
            self.send_fn(f"COLD_FLOW_{side}_{run_time_ms}")
            self.cf_status.setText(f"Running {side} for {run_time}s...")

    def _start_static_fire(self):
        run_time = self.sf_time.value()
        reply = QMessageBox.question(
            self, "Confirm Static Fire",
            f"Run static fire for {run_time}s?\n\nThis will ignite the engine.",
            QMessageBox.Yes | QMessageBox.No
        )
        if reply == QMessageBox.Yes:
            run_time_ms = run_time * 1000
            self.send_fn(f"STATIC_FIRE_{run_time_ms}")
            self.sf_status.setText(f"Running for {run_time}s...")

    def _toggle_logging(self):
        self.logging_active = not self.logging_active
        if self.logging_active:
            # This means logging just start.
            self.log_rows.clear()
            self.log_start_time = datetime.now()
            self.log_folder_path = self.parent_folder / self.log_start_time.strftime('%Y-%m-%d_%H-%M-%S')
            self.log_btn.setText("⏹  Stop Logging")
            self.log_btn.setStyleSheet("background:#da3633; color:#fff; border-radius:4px;")
        else:
            self.log_btn.setText("⏺  Start Logging")
            self.log_btn.setStyleSheet("background:#388bfd; color:#fff; border-radius:4px;")

    def _export_csv(self):
        if not self.log_rows:
            QMessageBox.information(self, "No Data", "No logged data to export.")
            return
        
        if self.log_folder_path is not None and self.log_start_time is not None:
            fname = self.log_folder_path / f"gse_log_{self.log_start_time.strftime('%Y%m%d_%H%M%S')}.csv"
        else:
            fname, _ = QFileDialog.getSaveFileName(
                self, "Save CSV", f"gse_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                "CSV Files (*.csv)"
            )
        
        if fname:
            with open(fname, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=self.log_rows[0].keys())
                writer.writeheader()
                writer.writerows(self.log_rows)
            QMessageBox.information(self, "Exported", f"Saved {len(self.log_rows)} rows to:\n{fname}")

    def _load_backup(self):
        """Utility to convert a raw pickle file to CSV (for backup data storage)."""
        if self.logging_active:
            QMessageBox.warning(self, "Logging Active", "Stop logging before loading backup data.")
            return
        
        if self.serial_thread and self.serial_thread.isRunning():
            QMessageBox.warning(self, "Connected", "Disconnect from serial port before loading backup data.")
            return
        
        folder_dir = QFileDialog.getExistingDirectory(self, "Open Backup Folder", "")
        if not folder_dir:
            return
        
        
        file_pre_num_string = "data-"
        file_post_num_string = ".pkl"
        file_pre_num_count = len(file_pre_num_string)
        file_post_num_count = len(file_post_num_string)

        file_names = os.listdir(folder_dir)
        numbers = [int(file[file_pre_num_count:-file_post_num_count])\
                   for file in file_names if \
                    file.endswith(file_post_num_string) and\
                    file.startswith(file_pre_num_string)]
    
        num_files = len(numbers)
        if num_files == 0:
            # TODO: Check this works and doesn't mess stuff up
            QMessageBox.warning(self, "No Data Files", "No valid backup files found in the selected folder.")
            return
        
        numbers = sorted(numbers)

        sorted_packets = [f"{file_pre_num_string}{num}{file_post_num_string}"\
                        for num in numbers]
        
        
        data_rows = []

        for file in sorted_packets:
            if (int(file[file_pre_num_count:-file_post_num_count]) % 100 == 0):
                print(f"Reading file: {file}/{num_files} ({(int(file[file_pre_num_count:-file_post_num_count])/num_files)*100:.2f}%)")
            with open(f"{folder_dir}/{file}", 'rb') as f:
                try:
                    raw_read = pickle.load(f)
                    if isinstance(raw_read, list):
                        data_rows.extend(raw_read)
                    else:
                        data_rows.append(raw_read)
                except Exception as e:
                    print(f"Error loading {file}: {e}")
                    continue

        start_time = folder_dir.split("/")[-1] # get filename from path
        self.log_rows.clear()
        self.log_start_time = datetime.strptime(start_time, '%Y-%m-%d_%H-%M-%S')
        self.log_rows = data_rows
        self.log_folder_path = Path(folder_dir)
        self.log_count_lbl.setText(f"{len(self.log_rows)} rows loaded")

        # self._backup_connect()
        if should_load_to_charts := QMessageBox.question(
            self, "Load to Charts?", "Data was loaded internally, and can be exported with the export CSV function.\
                Load backup data into live charts? (May be slow for large datasets)",
            QMessageBox.Yes | QMessageBox.No
        ) == QMessageBox.Yes:
            fake_time = 0
            for row in self.log_rows:
                new_row = {}
                for key, value in row.items():
                    if key == 'time_s':
                        new_row['time_s'] = fake_time
                        fake_time += 1
                    if key == 'ET_UP':
                        new_row['PT_EU'] = row['ET_UP']
                    if key == 'ET_DN':
                        new_row['PT_ED'] = row['ET_DN']
                    if key == 'NIT_UP':
                        new_row['PT_NU'] = row['NIT_UP']
                    if key == 'NIT_DN':
                        new_row['PT_ND'] = row['NIT_DN']
                    if key == 'LC1':
                        new_row['LC_E1'] = row['LC1']
                    if key == 'LC2':
                        new_row['LC_E2'] = row['LC2']
                    if key == 'LC3':
                        new_row['LC_E3'] = row['LC3']
                    if key == 'ET_TOTAL':
                        new_row['ET_TOTAL'] = row['ET_TOTAL']
                    if key == 'NLC1':
                        new_row['LC_N1'] = row['NLC1']
                    if key == 'NLC2':
                        new_row['LC_N2'] = row['NLC2']
                    if key == 'NLC3':
                        new_row['LC_N3'] = row['NLC3']
                    if key == 'NIT_TOTAL':
                        new_row['NIT_TOTAL'] = row['NIT_TOTAL']
                    if key == 'LC_T':
                        new_row['LC_T'] = row['LC_T']
                    if key == 'ERV':
                        new_row['ERV'] = row['ERV']
                    if key == 'EV':
                        new_row['EV'] = row['EV']
                    if key == 'NRV':
                        new_row['NRV'] = row['NRV']
                    if key == 'NV':
                        new_row['NV'] = row['NV']
                    try:
                        new_row[key] = float(value)
                    except ValueError:
                        new_row[key] = value
            
                self._on_data(new_row)

        QMessageBox.information(self, "Backup Loaded", f"Loaded {len(self.log_rows)} rows from:\n{folder_dir}")
        

    # ── Theme ─────────────────────────────────────────────────────
    def _apply_dark_theme(self):
        self.setStyleSheet("""
            QMainWindow, QWidget  { background:#0d1117; color:#c9d1d9; }
            QGroupBox             { border:1px solid #30363d; border-radius:6px;
                                    margin-top:10px; padding-top:6px; color:#8b949e; }
            QGroupBox::title      { subcontrol-origin:margin; left:8px; }
            QComboBox, QSpinBox   { background:#161b22; color:#c9d1d9;
                                    border:1px solid #30363d; border-radius:4px; padding:4px; }
            QPushButton           { background:#21262d; color:#c9d1d9;
                                    border:1px solid #30363d; border-radius:4px; padding:6px; }
            QPushButton:hover     { background:#30363d; }
            QLabel                { color:#c9d1d9; }
        """)

    def closeEvent(self, event):
        self._disconnect()
        event.accept()


# ─────────────────────────────────────────
#  Entry Point
# ─────────────────────────────────────────
if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setFont(QFont("Courier New", 9))
    win = GroundStation()
    win.show()
    sys.exit(app.exec_())