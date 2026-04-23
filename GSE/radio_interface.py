#!/usr/bin/env python3

import queue
import threading
import time
import tkinter as tk
from tkinter import font as tkfont
import serial

PORT = "COM9"          # change if needed
BAUD_RATE = 57600

ser = None
serial_thread = None
stop_event = threading.Event()
ui_queue = queue.Queue()


def now():
    return time.strftime("%H:%M:%S")


# ===================== SERIAL =====================

def serial_worker():
    global ser
    try:
        ser = serial.Serial(PORT, BAUD_RATE, timeout=0.1)
        time.sleep(2.0)

        ui_queue.put(("status", f"Connected to {PORT}", "#00aa66"))
        ui_queue.put(("log", f"[{now()}] --- serial port opened ---"))
        ui_queue.put(("connected", None))

        partial = b""

        while not stop_event.is_set():
            if ser.in_waiting:
                data = ser.read(ser.in_waiting)
                partial += data

                while b"\n" in partial:
                    line, partial = partial.split(b"\n", 1)
                    line = line.decode(errors="replace").rstrip("\r")
                    ui_queue.put(("log", f"[{now()}] {line}"))
            else:
                time.sleep(0.02)

    except Exception as e:
        ui_queue.put(("status", f"Error: {e}", "#cc3333"))

    finally:
        if ser and ser.is_open:
            ser.close()
        ui_queue.put(("disconnected", None))


def connect():
    global serial_thread
    stop_event.clear()
    btn_connect.config(state="disabled")
    serial_thread = threading.Thread(target=serial_worker, daemon=True)
    serial_thread.start()


def disconnect():
    stop_event.set()


def send(cmd):
    if ser and ser.is_open:
        ser.write((cmd + "\n").encode())
        ser.flush()
        log_local(f">>> {cmd}")


# ===================== UI =====================

def set_status(msg, color="#ccc"):
    status_label.config(text=msg, fg=color)


def log_local(msg):
    log_text.config(state="normal")
    log_text.insert("end", f"[{now()}] {msg}\n")
    log_text.see("end")
    log_text.config(state="disabled")


def poll_ui():
    while True:
        try:
            kind, *data = ui_queue.get_nowait()
        except queue.Empty:
            break

        if kind == "status":
            set_status(*data)

        elif kind == "log":
            log_text.config(state="normal")
            log_text.insert("end", data[0] + "\n")
            log_text.see("end")
            log_text.config(state="disabled")

        elif kind == "connected":
            enable_controls(True)
            disconnect_btn.config(state="normal")

        elif kind == "disconnected":
            enable_controls(False)
            disconnect_btn.config(state="disabled")
            btn_connect.config(state="normal")

    root.after(50, poll_ui)


def enable_controls(state):
    s = "normal" if state else "disabled"
    for b in all_buttons:
        b.config(state=s)


# ===================== COMMANDS =====================

COMMANDS = [
    ("ETHANE_VENT_ON", "Ethane Vent ON"),
    ("ETHANE_VENT_OFF", "Ethane Vent OFF"),
    ("NITROUS_VENT_ON", "Nitrous Vent ON"),
    ("NITROUS_VENT_OFF", "Nitrous Vent OFF"),

    ("ETHANE_RUN_ON", "Ethane Run ON"),
    ("ETHANE_RUN_OFF", "Ethane Run OFF"),
    ("NITROUS_RUN_ON", "Nitrous Run ON"),
    ("NITROUS_RUN_OFF", "Nitrous Run OFF"),

    ("ETHANE_A", "+90°"),
    ("ETHANE_S", "+10°"),
    ("ETHANE_D", "360°"),
    ("ETHANE_R", "Reset"),

    ("NITROUS_A", "+90°"),
    ("NITROUS_S", "+10°"),
    ("NITROUS_D", "360°"),
    ("NITROUS_R", "Reset"),

    ("START", "Start"),
    ("STOP", "Stop"),
    ("HELP", "Help"),
]


# ===================== BUILD UI =====================

root = tk.Tk()
root.title("RFD900x Control")
root.geometry("900x800")
root.configure(bg="#0a0f1a")

title_font = tkfont.Font(family="Courier New", size=16, weight="bold")
mono = tkfont.Font(family="Courier New", size=10)

# Header
tk.Label(root, text="RFD900x CONTROL",
         font=title_font, bg="#0a0f1a", fg="#00ffaa").pack(pady=8)

topbar = tk.Frame(root, bg="#0a0f1a")
topbar.pack(fill="x", padx=10)

btn_connect = tk.Button(topbar, text="CONNECT",
                        command=connect, bg="#1f2937", fg="white")
btn_connect.pack(side="left", padx=5)

disconnect_btn = tk.Button(topbar, text="DISCONNECT",
                           command=disconnect, state="disabled",
                           bg="#552233", fg="white")
disconnect_btn.pack(side="left", padx=5)

status_label = tk.Label(topbar, text="Not connected",
                        bg="#0a0f1a", fg="#888")
status_label.pack(side="right")

# Controls (compact)
control_frame = tk.Frame(root, bg="#0a0f1a")
control_frame.pack(fill="x", padx=10, pady=5)

all_buttons = []

for i, (cmd, label) in enumerate(COMMANDS):
    b = tk.Button(control_frame,
                  text=label,
                  command=lambda c=cmd: send(c),
                  width=14,
                  bg="#1f2937",
                  fg="white",
                  state="disabled")
    b.grid(row=i // 6, column=i % 6, padx=4, pady=4)
    all_buttons.append(b)

# BIG LOG AREA
log_frame = tk.Frame(root, bg="#0a0f1a")
log_frame.pack(fill="both", expand=True, padx=10, pady=10)

scrollbar = tk.Scrollbar(log_frame)
scrollbar.pack(side="right", fill="y")

log_text = tk.Text(
    log_frame,
    font=mono,
    bg="#050a12",
    fg="#00ffaa",
    insertbackground="white",
    yscrollcommand=scrollbar.set
)
log_text.pack(fill="both", expand=True)

scrollbar.config(command=log_text.yview)

# Start loop
enable_controls(False)
root.after(50, poll_ui)
root.mainloop()