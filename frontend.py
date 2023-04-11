from tkinter import *
from tkinter.filedialog import askopenfile
import tkinter as tk
import customtkinter as ctk
import os

ctk.set_appearance_mode("System")
ctk.set_default_color_theme("green")

class App(ctk.CTk):
    def __init__(self):
        super().__init__()
#instantiate the main app. this will be cleared once the button is pressed.
        self.title("DeGenPrime")
        self.geometry("500x300")
        frame = Frame(self)
        frame.grid()

        #event handling for when the button is pressed.
        def handle_click(event):
            
            # configure grid layout (4x4)
            self.grid_columnconfigure((0,1,2,3), weight=1)
            self.grid_rowconfigure((0, 1, 2), weight=1)

            #delete old info
            self.btn_greeting.destroy()
            self.lbl_greeting.destroy()

            def find_input_file():
                file = askopenfile(mode ='r', filetypes =[('Acceptable Filetypes', '*.fasta,*.fna,*.clust,*.faa')])
                if file:
                    inputpath = os.path.abspath(file.name)
                    self.lbl_inputpath = ctk.CTkLabel(self.sidebar_frame, text="Your selected file is: " + str(inputpath),wraplength = 140,justify="center")
                    self.lbl_inputpath.grid(row=2, column=0, padx=20)
                    filecontent = file.read()

            #create the frame for the three side buttons -- browse, run, and save
            self.sidebar_frame = ctk.CTkFrame(self, width=140, corner_radius=0)
            self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
            self.sidebar_frame.grid_rowconfigure(4, weight=1)
            self.btn_browse = ctk.CTkButton(self.sidebar_frame, text = "BROWSE FILES", command = find_input_file)
            self.btn_browse.grid(row=0, column=0, padx=20, pady=30)
            self.btn_run = ctk.CTkButton(self.sidebar_frame, text = "RUN PROGRAM")
            self.btn_run.grid(row=1, column=0, padx=20, pady=30)
            self.btn_save = ctk.CTkButton(self.sidebar_frame, text = "SAVE SETTINGS")
            self.btn_save.grid(row=2, column=0, padx=20, pady=30)

            self.grid_rowconfigure(0, weight=1)
            self.grid_columnconfigure((0, 1, 2), weight=1)

            #create the frame for tabview (settings and advanced settings)
            self.settings_tabview = ctk.CTkTabview(self, width=250)
            self.settings_tabview.grid(row=0, column=1, padx=(20, 0), pady=(20, 0), sticky="nsew")
            self.settings_tabview.add("Settings")
            self.settings_tabview.add("Advanced Settings")

            self.amplicon = ctk.CTkEntry(self.settings_tabview.tab("Settings"), placeholder_text="Amplicon Length")
            self.amplicon.grid(row = 0, column = 0,padx=20, pady=(20, 10))
            self.ampl_label = ctk.CTkLabel(self.settings_tabview.tab("Settings"), text= "Base Pairs")
            self.ampl_label.grid(row = 0, column = 2,pady=(20, 10))
            self.DNAconc = ctk.CTkEntry(self.settings_tabview.tab("Settings"), placeholder_text="DNA Concentration")
            self.DNAconc.grid(row = 1, column = 0,padx=20, pady=(20, 10))
            self.DNAconc_label = ctk.CTkLabel(self.settings_tabview.tab("Settings"), text= "mMol")
            self.DNAconc_label.grid(row = 1, column = 2,pady=(20, 10))
            self.temp = ctk.CTkEntry(self.settings_tabview.tab("Settings"), placeholder_text="Temperature")
            self.temp.grid(row = 2, column = 0,padx=20, pady=(20, 10))
            self.temp_label = ctk.CTkLabel(self.settings_tabview.tab("Settings"), text= "C°")
            self.temp_label.grid(row = 2, column = 2,pady=(20, 10))

            self.primer_conc = ctk.CTkEntry(self.settings_tabview.tab("Advanced Settings"), placeholder_text="Primer Concentration")
            self.primer_conc.grid(row = 0, column = 0,padx=20, pady=(20, 10))
            self.primer_conc_label = ctk.CTkLabel(self.settings_tabview.tab("Advanced Settings"), text= "nMol")
            self.primer_conc_label.grid(row = 0, column = 2,pady=(20, 10))
            self.mono_conc = ctk.CTkEntry(self.settings_tabview.tab("Advanced Settings"), placeholder_text="Monovalent Ion Concentration")
            self.mono_conc.grid(row = 1, column = 0,padx=20, pady=(20, 10))
            self.mono_conc_label = ctk.CTkLabel(self.settings_tabview.tab("Advanced Settings"), text= "mMol")
            self.mono_conc_label.grid(row = 1, column = 2,pady=(20, 10))

        self.sidebar_frame = ctk.CTkFrame(self, width = 500)
        self.lbl_greeting = ctk.CTkLabel(self,text="Welcome to DeGenPrime. Our goal is to maximize your primer design experience while minimizing your headaches. "
                                         + "Input either an unaligned or a prealigned .fasta or .fna file, and allow DeGenPrime to recommend the best primers for your PCR use."
                                          ,wraplength = 450, justify= "center")
        self.lbl_greeting.grid(padx=30, pady=50)
        self.btn_greeting = ctk.CTkButton(self,text="ENTER")
        self.btn_greeting.bind("<Button-1>", handle_click)
        self.btn_greeting.grid(padx=30, pady=10)

if __name__ == "__main__":
    app = App()
    app.mainloop()