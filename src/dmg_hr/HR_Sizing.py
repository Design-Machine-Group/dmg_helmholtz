import tkinter as tk
from tkinter import *
import numpy as np
import math
import os

global entries, fields, c
entries = []
fields = ['Ln - Neck Length (cm): ', 'Rn - Neck Radius (cm): ', 'Rn\' - Neck Radius (cm): ', 'Rc - Cavity Radius (cm): ', 'Lc - Cavity Depth (cm): ']
c = 343 # Speed of sound in air


def eval_freq():
    global radio_val, res_freq, entries, c
    ln = float(entries[0][1].get()) * 0.01
    rn1 = float(entries[1][1].get()) * 0.01
    rn2 = float(entries[2][1].get()) * 0.01
    rc = float(entries[3][1].get()) * 0.01
    lc = float(entries[4][1].get()) * 0.01

    if len([x for x in [ln,rn1,rc,lc] if x > 0]) == 4:
        # Length end corrections
        ln = ln + 1.7*rn1

        # Surface area of the opening
        if radio_val.get() == 4:
            S = np.pi * rn1 * rn2
        else:
            S = np.pi * rn1 * rn1

        # Volume of the cavity
        if radio_val.get() == 1:
            V = np.pi * rc * rc * lc
        elif radio_val.get() == 2:
            V = rc*rc*lc
        elif radio_val.get() == 3:
            V = (np.pi * rc * rc * lc) - (np.pi * rn1 * rn1 * ln)
        elif radio_val.get() == 4:
            V = np.pi * rc * rc * lc

        freq = (c/(2*np.pi)) * math.sqrt(S/(V*ln))
        res_freq.delete(0, END)
        res_freq.insert(0, freq)

        return True

def makeform(root):
    global entries, fields
    i = 1
    for field in fields:
        row = tk.Frame(root)
        lab = tk.Label(row, width=25, text=field)

        sv = StringVar()
        ent = tk.Entry(row, textvariable=sv, validate="focusout", validatecommand=eval_freq)

        ent.insert(0, "0")
        row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        lab.pack(side=tk.LEFT)
        ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
        entries.append((lab, ent))
        i+=1
    return entries

def show_fields_circular():
    global fields
    entries[0][0]['text'] = fields[0]
    entries[1][0]['text'] = fields[1]
    entries[2][1].config(state='disabled')
    entries[3][0]['text'] = fields[3]
    entries[4][0]['text'] = fields[4]
    eval_freq()

def show_fields_rectangular():
    global fields
    entries[0][0]['text'] = fields[0]
    entries[1][0]['text'] = fields[1]
    entries[2][1].config(state='disabled')
    entries[3][0]['text'] = 'Wc - Cavity Width (cm): '
    entries[4][0]['text'] = fields[4]
    eval_freq()

def show_fields_embedded():
    global fields
    entries[0][0]['text'] = ' Ln - Embedded Length (cm): '
    entries[1][0]['text'] = fields[1]
    entries[2][1].config(state='disabled')
    entries[3][0]['text'] = fields[3]
    entries[4][0]['text'] = fields[4]
    eval_freq()

def show_fields_conical():
    global fields
    entries[0][0]['text'] = fields[0]
    entries[1][0]['text'] = 'Rn - Outer Neck Radius (cm): '
    entries[2][1].config(state='normal')
    entries[2][0]['text'] = 'Rn\' - Inner Neck Radius (cm): '
    entries[3][0]['text'] = fields[3]
    entries[4][0]['text'] = fields[4]
    eval_freq()

def generate_graph():
    global radio_val, entries, c





if __name__ == '__main__':
    global radio_val, res_freq
    master = Tk()
    master.title("Helmholtz resontor design")
    master.resizable()

    # Basic layout
    ents = makeform(master)
    entries[2][1].config(state='disabled')

    # Adding buttons
    row1 = tk.Frame(master)
    radio_val = tk.IntVar()
    radio_val.set(1)
    radio1 = tk.Radiobutton(row1, text='Circular', command=lambda: show_fields_circular(), variable=radio_val, value=1)
    radio1.pack(side=tk.LEFT, padx=5, pady=15)
    radio2 = tk.Radiobutton(row1, text='Rectangular', command=lambda: show_fields_rectangular(), variable=radio_val, value=2)
    radio2.pack(side=tk.LEFT, padx=5, pady=15)
    radio3 = tk.Radiobutton(row1, text='Embedded neck', command=lambda: show_fields_embedded(), variable=radio_val, value=3)
    radio3.pack(side=tk.LEFT, padx=5, pady=15)
    radio4 = tk.Radiobutton(row1, text='Conical neck', command=lambda: show_fields_conical(), variable=radio_val, value=4)
    radio4.pack(side=tk.LEFT, padx=5, pady=15)
    row1.pack()

    # Adding the resonant freq text box
    row2 = tk.Frame(master)
    w = tk.Label(row2, width=25, text='Resonant Frequency: ')
    w.pack(side=tk.LEFT)

    res_freq = Entry(row2)
    res_freq.insert(0, "0")
    res_freq.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
    row2.pack(fill=tk.X, padx=15, pady=15)

    # Adding the Graph button
    Button(master, text="Generate Graph", command=generate_graph).pack(padx=15, pady=25)

    mainloop()