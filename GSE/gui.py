import sys
import csv
import threading
from collections import deque
from datetime import datetime

import serial
import serial.tools.list_ports


from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QComboBox, QSpinBox,
    QGroupBox, QSizePolicy, QFileDialog, QMessageBox, QFrame, QTextEdit
)
from PyQt5.QtCore import Qt, QTimer, pyqtSignal, QObject, QThread
from PyQt5.QtGui import QFont, QColor, QPalette

import pyqtgraph as pg


# ─────────────────────────────────────────
#  Serial Worker (runs in background thread)
# ─────────────────────────────────────────
class SerialWorker(QObject):
    data_received = pyqtSignal(dict)
    raw_received = pyqtSignal(str)  # Raw serial line for monitor
    connection_lost = pyqtSignal()

    def __init__(self, port, baud=9600):
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

        while self._running:
            try:
                raw = self.ser.readline().decode('utf-8', errors='ignore').strip()
                if raw.startswith("DATA|"):
                    parsed = self._parse(raw)
                    if parsed:
                        self.data_received.emit(parsed)
                    self.raw_received.emit(raw)  # Emit raw line for monitor
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
        """Parse DATA|millis|KEY:VAL|KEY:VAL|... into a dict."""
        try:
            parts = line.split('|')
            # parts[0] = "DATA", parts[1] = millis
            result = {'millis': int(parts[1])}
            for part in parts[2:]:
                if ':' in part:
                    k, v = part.split(':', 1)
                    result[k] = float(v)
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
        self.thread = None
        self.send_fn = lambda cmd: None  # no-op until connected

        # Data history for charts
        self.t_hist    = deque(maxlen=self.HISTORY_LEN)
        self.et_up_hist = deque(maxlen=self.HISTORY_LEN)
        self.et_dn_hist = deque(maxlen=self.HISTORY_LEN)
        self.nit_up_hist = deque(maxlen=self.HISTORY_LEN)
        self.nit_dn_hist = deque(maxlen=self.HISTORY_LEN)
        self.lc_total_hist = deque(maxlen=self.HISTORY_LEN)

        # CSV logging
        self.log_rows = []
        self.logging_active = False

        # Serial monitor
        self.serial_monitor_lines = deque(maxlen=500)

        self._build_ui()

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
        left.addStretch()

        right = QVBoxLayout()
        right.addWidget(self._build_valve_group())
        right.addWidget(self._build_cold_flow_group())
        right.addWidget(self._build_export_group())
        right.addWidget(self._build_serial_monitor_group())
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

        self.lc_e1 = SensorLabel("Ethane LC1",  "kg")
        self.lc_e2 = SensorLabel("Ethane LC2",  "kg")
        self.lc_e3 = SensorLabel("Ethane LC3",  "kg")
        self.lc_et = SensorLabel("Ethane TOTAL","kg")
        self.lc_n1 = SensorLabel("Nitrous LC1", "kg")
        self.lc_n2 = SensorLabel("Nitrous LC2", "kg")
        self.lc_n3 = SensorLabel("Nitrous LC3", "kg")
        self.lc_nt = SensorLabel("Nitrous TOTAL","kg")

        for i, w in enumerate([self.lc_e1, self.lc_e2, self.lc_e3, self.lc_et]):
            grid.addWidget(w, 0, i)
        for i, w in enumerate([self.lc_n1, self.lc_n2, self.lc_n3, self.lc_nt]):
            grid.addWidget(w, 1, i)
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

        self.log_count_lbl = QLabel("0 rows logged")
        self.log_count_lbl.setFont(QFont("Courier New", 8))
        self.log_count_lbl.setStyleSheet("color:#8b949e;")

        layout.addWidget(self.log_btn)
        layout.addWidget(export_btn)
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

        clear_btn = QPushButton("Clear")
        clear_btn.setFixedHeight(28)
        clear_btn.setFont(QFont("Courier New", 8))
        clear_btn.clicked.connect(self.serial_monitor.clear)

        layout.addWidget(self.serial_monitor)
        layout.addWidget(clear_btn)
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

        self.lc_chart = pg.PlotWidget(title="Propellant Mass (kg)")
        self.lc_chart.addLegend()
        self.lc_chart.showGrid(x=True, y=True, alpha=0.2)
        self.curve_lc_et  = self.lc_chart.plot(pen=pg.mkPen('#3fb950', width=2), name="Ethane")
        self.curve_lc_nit = self.lc_chart.plot(pen=pg.mkPen('#d29922', width=2), name="Nitrous")

        layout.addWidget(self.pt_chart)
        layout.addWidget(self.lc_chart)
        return grp

    # ── Connection ────────────────────────────────────────────────
    def _refresh_ports(self):
        self.port_combo.clear()
        ports = [p.device for p in serial.tools.list_ports.comports()]
        self.port_combo.addItems(ports if ports else ["No ports found"])

    def _toggle_connection(self):
        if self.thread and self.thread.isRunning():
            self._disconnect()
        else:
            self._connect()

    def _connect(self):
        port = self.port_combo.currentText()
        if not port or port == "No ports found":
            QMessageBox.warning(self, "No Port", "Select a valid serial port.")
            return

        self.worker = SerialWorker(port)
        self.thread = QThread()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.start)
        self.worker.data_received.connect(self._on_data)
        self.worker.raw_received.connect(self._on_raw_data)
        self.worker.connection_lost.connect(self._disconnect)
        self.send_fn = self.worker.send
        self.thread.start()

        self.connect_btn.setText("Disconnect")
        self.connect_btn.setStyleSheet("background:#da3633; color:#fff; border-radius:4px;")
        self.status_lbl.setText("● Connected")
        self.status_lbl.setStyleSheet("color:#3fb950;")

    def _disconnect(self):
        if self.worker:
            self.worker.stop()
        if self.thread:
            self.thread.quit()
            self.thread.wait()
        self.send_fn = lambda cmd: None
        self.connect_btn.setText("Connect")
        self.connect_btn.setStyleSheet("background:#238636; color:#fff; border-radius:4px;")
        self.status_lbl.setText("● Disconnected")
        self.status_lbl.setStyleSheet("color:#ff4466;")
        self.serial_monitor.clear()

    # ── Data Handler ──────────────────────────────────────────────
    def _on_data(self, state: dict):
        t = state.get('millis', 0) / 1000.0

        # PT readouts
        et_up  = state.get('PT_EU',  0.0)
        et_dn  = state.get('PT_ED',  0.0)
        nit_up = state.get('PT_NU', 0.0)
        nit_dn = state.get('PT_ND', 0.0)
        self.pt_et_up.update_value(et_up)
        self.pt_et_dn.update_value(et_dn)
        self.pt_nit_up.update_value(nit_up)
        self.pt_nit_dn.update_value(nit_dn)

        # LC readouts
        lc1 = state.get('LC_E1', 0.0)
        lc2 = state.get('LC_E2', 0.0)
        lc3 = state.get('LC_E3', 0.0)
        et_total = lc1 + lc2 + lc3
        n1  = state.get('LC_N1', 0.0)
        n2  = state.get('LC_N2', 0.0)
        n3  = state.get('LC_N3', 0.0)
        nit_total = n1 + n2 + n3
        self.lc_e1.update_value(lc1)
        self.lc_e2.update_value(lc2)
        self.lc_e3.update_value(lc3)
        self.lc_et.update_value(et_total)
        self.lc_n1.update_value(n1)
        self.lc_n2.update_value(n2)
        self.lc_n3.update_value(n3)
        self.lc_nt.update_value(nit_total)

        # Valve states from firmware flags
        self.btn_erv.set_state(bool(state.get('ERV', 0)))
        self.btn_ev.set_state( bool(state.get('EV',  0)))
        self.btn_nrv.set_state(bool(state.get('NRV', 0)))
        self.btn_nv.set_state( bool(state.get('NV',  0)))

        # Charts
        self.t_hist.append(t)
        self.et_up_hist.append(et_up)
        self.et_dn_hist.append(et_dn)
        self.nit_up_hist.append(nit_up)
        self.nit_dn_hist.append(nit_dn)
        self.lc_total_hist.append(et_total)

        tl = list(self.t_hist)
        self.curve_et_up.setData(tl, list(self.et_up_hist))
        self.curve_et_dn.setData(tl, list(self.et_dn_hist))
        self.curve_nit_up.setData(tl, list(self.nit_up_hist))
        self.curve_nit_dn.setData(tl, list(self.nit_dn_hist))
        self.curve_lc_et.setData(tl, list(self.lc_total_hist))

        # CSV logging
        if self.logging_active:
            self.log_rows.append({
                'time_s': t, 'ET_UP': et_up, 'ET_DN': et_dn,
                'NIT_UP': nit_up, 'NIT_DN': nit_dn,
                'LC1': lc1, 'LC2': lc2, 'LC3': lc3, 'ET_TOTAL': et_total,
                'NLC1': n1, 'NLC2': n2, 'NLC3': n3, 'NIT_TOTAL': nit_total,
                'ERV': state.get('ERV',0), 'EV': state.get('EV',0),
                'NRV': state.get('NRV',0), 'NV': state.get('NV',0),
            })
            self.log_count_lbl.setText(f"{len(self.log_rows)} rows logged")

    def _on_raw_data(self, raw: str):
        """Append raw serial line to monitor."""
        self.serial_monitor.append(raw)
        # Auto-scroll to bottom
        self.serial_monitor.verticalScrollBar().setValue(
            self.serial_monitor.verticalScrollBar().maximum()
        )

    # ── Controls ──────────────────────────────────────────────────
    def _emergency_stop(self):
        for cmd in ["ETHANE_RUN_OFF", "NITROUS_RUN_OFF",
                    "ETHANE_VENT_ON", "NITROUS_VENT_ON"]:
            self.send_fn(cmd)

    def _start_cold_flow(self):
        side     = self.cf_side.currentText()
        run_time = self.cf_time.value()
        reply = QMessageBox.question(
            self, "Confirm Cold Flow",
            f"Run cold flow on {side} side for {run_time}s?\n\nThis will open valves.",
            QMessageBox.Yes | QMessageBox.No
        )
        if reply == QMessageBox.Yes:
            self.send_fn(f"COLD_FLOW_{side}_{run_time}")
            self.cf_status.setText(f"Running {side} for {run_time}s...")

    def _toggle_logging(self):
        self.logging_active = not self.logging_active
        if self.logging_active:
            self.log_rows.clear()
            self.log_btn.setText("⏹  Stop Logging")
            self.log_btn.setStyleSheet("background:#da3633; color:#fff; border-radius:4px;")
        else:
            self.log_btn.setText("⏺  Start Logging")
            self.log_btn.setStyleSheet("background:#388bfd; color:#fff; border-radius:4px;")

    def _export_csv(self):
        if not self.log_rows:
            QMessageBox.information(self, "No Data", "No logged data to export.")
            return
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