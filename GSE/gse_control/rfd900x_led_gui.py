#!/usr/bin/env python3
"""
rfd900x_led_gui.py
Tkinter GUI to control Arduino LED over RFD900x radio on COM3.

Usage:
    pip install pyserial
    python rfd900x_led_gui.py
"""

import serial
import threading
import time
import tkinter as tk
from tkinter import font as tkfont

PORT      = "COM3"
BAUD_RATE = 57600

# ── Serial state ───────────────────────────────────────────────────────────────
ser = None
led_state = False


def connect():
    global ser
    try:
        ser = serial.Serial(PORT, BAUD_RATE, timeout=1)
        time.sleep(2)
        ser.reset_input_buffer()
        set_status("Connected", "#00ff99")
        btn_on.config(state="normal")
        btn_off.config(state="normal")
        btn_connect.config(state="disabled")
        threading.Thread(target=reader, daemon=True).start()
    except serial.SerialException as e:
        set_status(f"Error: {e}", "#ff4444")


def reader():
    """Background thread — reads replies from Arduino."""
    while ser and ser.is_open:
        try:
            line = ser.readline().decode(errors="replace").strip()
            if line:
                log(f"← {line}")
                if "on" in line.lower():
                    set_led_indicator(True)
                elif "off" in line.lower():
                    set_led_indicator(False)
        except Exception:
            break


def send(cmd):
    if ser and ser.is_open:
        ser.write((cmd + "\n").encode())
        ser.flush()
        log(f"→ {cmd}")
    else:
        set_status("Not connected", "#ff4444")


def led_on():
    send("1")

def led_off():
    send("0")


def set_led_indicator(state: bool):
    color = "#00ff99" if state else "#1a1a2e"
    glow  = "#00ff9966" if state else "#00000000"
    label = "ON" if state else "OFF"
    root.after(0, lambda: led_circle.config(bg=color))
    root.after(0, lambda: led_label.config(text=f"LED  {label}",
                                            fg="#00ff99" if state else "#444466"))


def set_status(msg, color="#aaaacc"):
    root.after(0, lambda: status_label.config(text=msg, fg=color))


def log(msg):
    def _log():
        log_text.config(state="normal")
        log_text.insert("end", msg + "\n")
        log_text.see("end")
        log_text.config(state="disabled")
    root.after(0, _log)


# ── Build GUI ──────────────────────────────────────────────────────────────────
root = tk.Tk()
root.title("RFD900x LED Control")
root.configure(bg="#0d0d1a")
root.resizable(False, False)

W = 360
root.geometry(f"{W}x520")

mono  = tkfont.Font(family="Courier New", size=11)
title_font = tkfont.Font(family="Courier New", size=15, weight="bold")
small_font = tkfont.Font(family="Courier New", size=9)

# ── Title ──
tk.Label(root, text="RFD900x", font=title_font,
         bg="#0d0d1a", fg="#00ff99").pack(pady=(24, 0))
tk.Label(root, text="LED CONTROLLER", font=tkfont.Font(family="Courier New", size=9),
         bg="#0d0d1a", fg="#444466").pack()

# ── Divider ──
tk.Frame(root, bg="#1e1e3a", height=1).pack(fill="x", padx=24, pady=16)

# ── LED indicator ──
indicator_frame = tk.Frame(root, bg="#0d0d1a")
indicator_frame.pack(pady=8)

led_circle = tk.Label(indicator_frame, width=4, height=2,
                      bg="#1a1a2e", relief="flat", bd=0)
led_circle.pack(side="left", padx=(0, 12))

led_label = tk.Label(indicator_frame, text="LED  OFF",
                     font=tkfont.Font(family="Courier New", size=13, weight="bold"),
                     bg="#0d0d1a", fg="#444466")
led_label.pack(side="left")

# ── Buttons ──
btn_frame = tk.Frame(root, bg="#0d0d1a")
btn_frame.pack(pady=20)

btn_on = tk.Button(btn_frame, text="ON", width=8, height=2,
                   font=tkfont.Font(family="Courier New", size=13, weight="bold"),
                   bg="#003322", fg="#00ff99", activebackground="#005533",
                   activeforeground="#00ff99", relief="flat", bd=0,
                   cursor="hand2", state="disabled",
                   command=led_on)
btn_on.grid(row=0, column=0, padx=10)

btn_off = tk.Button(btn_frame, text="OFF", width=8, height=2,
                    font=tkfont.Font(family="Courier New", size=13, weight="bold"),
                    bg="#220011", fg="#ff4466", activebackground="#440022",
                    activeforeground="#ff4466", relief="flat", bd=0,
                    cursor="hand2", state="disabled",
                    command=led_off)
btn_off.grid(row=0, column=1, padx=10)

# ── Connect button ──
btn_connect = tk.Button(root, text=f"CONNECT  {PORT}", width=22, height=1,
                        font=tkfont.Font(family="Courier New", size=11),
                        bg="#1a1a2e", fg="#aaaacc", activebackground="#2a2a4e",
                        activeforeground="#ffffff", relief="flat", bd=0,
                        cursor="hand2", command=connect)
btn_connect.pack(pady=4)

# ── Status ──
status_label = tk.Label(root, text="Not connected", font=small_font,
                         bg="#0d0d1a", fg="#444466")
status_label.pack(pady=(4, 12))

# ── Divider ──
tk.Frame(root, bg="#1e1e3a", height=1).pack(fill="x", padx=24, pady=4)

# ── Log ──
tk.Label(root, text="LOG", font=small_font,
         bg="#0d0d1a", fg="#333355").pack(anchor="w", padx=24, pady=(8, 2))

log_frame = tk.Frame(root, bg="#0d0d1a")
log_frame.pack(fill="both", expand=True, padx=24, pady=(0, 24))

log_text = tk.Text(log_frame, height=7, font=mono,
                   bg="#080812", fg="#4444aa", insertbackground="#00ff99",
                   relief="flat", bd=0, state="disabled",
                   wrap="word")
log_text.pack(fill="both", expand=True)

root.mainloop()
