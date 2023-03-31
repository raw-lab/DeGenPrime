from tkinter import *
from tkinter.filedialog import askopenfile
import tkinter as tk
import customtkinter as ctk
import os

#instantiate the main window. this will be cleared once the button is pressed.
window = ctk.CTk()
window.title("DeGenPrime")
window.geometry("500x300")
frame = Frame(window)
frame.pack()


#event handling for when the button is pressed.
def handle_click(event):
    lbl_greeting.destroy()
    btn_greeting.destroy()
    window.geometry("700x500")
    btn_browsefiles.place(x=10,y=10)

#entry screen. this will describe to the user the purpose of DeGenPrime
lbl_greeting = ctk.CTkLabel(text="Welcome to DeGenPrime. Our goal is to maximize your primer design experience while minimizing your headaches. Input either an unaligned or a prealigned .fasta or .fna file, and allow DeGenPrime to recommend the best primers for your PCR use.", master = window,wraplength = 450, justify= "center")
lbl_greeting.pack(fill=BOTH, expand=True)
btn_greeting = ctk.CTkButton(master = window,text="ENTER")
btn_greeting.bind("<Button-1>", handle_click)
btn_greeting.pack()

#event handling for the user to input their file
def find_input_file():
    file = askopenfile(mode ='r', filetypes =[('Acceptable Filetypes', '*.fasta,*.fna,*.clust,*.faa')])
    if file:
        inputpath = os.path.abspath(file.name)
        lbl_inputpath = ctk.CTkLabel(master = window, text="Your selected file is: " + str(inputpath)).place(x=10,y=50)
        filecontent = file.read()

        #somewhere in this event handling scheme we will have to pass the file to C++. figure out how to do this.


btn_browsefiles = ctk.CTkButton(master=window,text="Browse Files", command=find_input_file)



window.mainloop()