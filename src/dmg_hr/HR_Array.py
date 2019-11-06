import tkinter as tk
from tkinter import *
import numpy as np
import math
import os

global entries
entries = {}

def makeform(root):
    global entries
    fields = ['Source position: ', 'Receiver position: ', 'Array position: ']
    for field in fields:
        row = tk.Frame(root)
        tk.Label(row, width=20, text=field).pack(side=tk.LEFT)
        val = Entry(row)
        val.insert(0, "0,0,0")
        val.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
        row.pack(fill=tk.X, padx=10, pady=15)
        entries.update({field: val})

    row = tk.Frame(root)
    fields = ['M resonators in X-axis', 'N resonators in Y-axis']
    for field in fields:
        tk.Label(row, width=20, text=field).pack(side=tk.LEFT)
        val = Entry(row)
        val.insert(0, "0")
        val.pack(side=tk.LEFT, padx=10, pady=1)
        entries.update({field: val})
    row.pack(fill=tk.X, padx=3, pady=3)

    row = tk.Frame(root)
    fields = ['Dist (necks - X) in cm: ', 'Dist (necks - Y) in cm: ']
    for field in fields:
        tk.Label(row, width=20, text=field).pack(side=tk.LEFT)
        val = Scale(row, from_=0, to=50, orient=HORIZONTAL)
        val.pack(side=tk.LEFT, padx=10, pady=1)
        entries.update({field: val})
    row.pack(fill=tk.X, padx=3, pady=3)

def generate_graph():
    global entries




if __name__ == '__main__':
    master = Tk()
    master.title("Helmholtz resontor design")
    master.resizable()

    # Basic layout
    makeform(master)

    # Adding the Graph button
    Button(master, text="Generate Graph", command=generate_graph).pack(padx=15, pady=25)

    mainloop()