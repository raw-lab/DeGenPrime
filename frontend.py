import tkinter as tk
from tkinter import *
from tkinter.filedialog import askopenfile
import os

#instantiate the main window. this will be cleared once the button is pressed.
window = tk.Tk()
window.title("DeGenPrime")
window.geometry("250x250")
frame = Frame(window)
frame.pack()


#event handling for when the button is pressed.
def handle_click(event):
    lbl_greeting.destroy()
    btn_greeting.destroy()
    window.geometry("500x500")
    btn_browsefiles.pack()

#entry screen. this will describe to the user the purpose of DeGenPrime
lbl_greeting = tk.Label(text="Welcome to DeGenPrime. Our goal is to maximize your primer design experience while minimizing your headaches. Input either an unaligned or a prealigned .fasta or .fna file, and allow DeGenPrime to recommend the best primers for your PCR use.")
lbl_greeting.pack(fill=BOTH, expand=True)
btn_greeting = tk.Button(text="ENTER")
btn_greeting.bind("<Button-1>", handle_click)
btn_greeting.pack()

#event handling for the user to input their file
def find_input_file():
    file = askopenfile(mode ='r', filetypes =[('Acceptable Filetypes', '*.fasta,*.fna,*.clust,*.faa')])
    if file:
        inputpath = os.path.abspath(file.name)
        lbl_inputpath = tk.Label(text="Your selected file is" + str(inputpath)).pack()

btn_browsefiles = tk.Button(text="Browse Files", command=find_input_file)



window.mainloop()