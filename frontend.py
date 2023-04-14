from tkinter import *
from tkinter.filedialog import askopenfile
import tkinter as tk
import customtkinter as ctk
import os, csv

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

            #event handling for finding the input file
            def find_input_file():
                file = askopenfile(mode ='r', filetypes =[('Acceptable Filetypes', '*.fasta *.fna *.clust *.faa')])
                if file:
                    inputpath = os.path.abspath(file.name)
                    self.lbl_inputpath = ctk.CTkLabel(self, text="Your selected file is: " + str(inputpath),wraplength = 300,justify="center")
                    self.lbl_inputpath.grid(row=1, column=1, padx = 10, pady=10)
                    filecontent = file.read()

            def save_settings():
                amplicon_get = self.amplicon.get()
                if(int(amplicon_get) < 0): amplicon_get = 0
                DNAconc_get = self.DNAconc.get()
                if(int(DNAconc_get) < 0): DNAconc_get = 50
                temp_get = self.temp.get()
                if(int(temp_get) < 50): temp_get = 50
                pc_get = self.pc.get()
                if(int(pc_get) < 0): pc_get = 50
                mc_get = self.mc.get()
                if(int(mc_get) < 0): mc_get = 50

                with open('SETTINGS.csv', 'w', newline='') as file:
                    writer = csv.writer(file)
                    field = ["Amplicon_Length", "DNA_Concentration", "Temperature", "Primer_Concentration", "Monovalent_Ion_Concentration"]
    
                    writer.writerow(field)
                    writer.writerow([amplicon_get, DNAconc_get, temp_get, pc_get, mc_get])

            #create the frame for the three side buttons -- browse, run, and save
            self.sidebar_frame = ctk.CTkFrame(self, width=140, corner_radius=0)
            self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
            self.sidebar_frame.grid_rowconfigure(4, weight=1)
            self.btn_browse = ctk.CTkButton(self.sidebar_frame, text = "BROWSE FILES", command = find_input_file)
            self.btn_browse.grid(row=0, column=0, padx=20, pady=(40,30))
            self.btn_run = ctk.CTkButton(self.sidebar_frame, text = "RUN PROGRAM")
            self.btn_run.grid(row=2, column=0, padx=20, pady=30)
            self.btn_save = ctk.CTkButton(self.sidebar_frame, text = "SAVE SETTINGS", command = save_settings)
            self.btn_save.grid(row=3, column=0, padx=20, pady=30)

            self.grid_rowconfigure((0,1), weight=1)
            self.grid_columnconfigure((0, 1, 2), weight=1)

            #create the frame for tabview (settings and advanced settings)
            self.settings_tabview = ctk.CTkTabview(self, width=250)
            self.settings_tabview.grid(row=0, column=1, padx=(20, 0), pady=(1, 0), sticky="nsew")
            self.settings_tabview.add("Settings")
            self.settings_tabview.add("Advanced Settings")

            ### REGULAR SETTINGS ###
            self.ampl_label = ctk.CTkLabel(self.settings_tabview.tab("Settings"), text= "Minimum Amplicon Length")
            self.ampl_label.grid(row = 0, column = 0,pady=2.5)
            self.amplicon = ctk.CTkEntry(self.settings_tabview.tab("Settings"), placeholder_text="Min. Amplicon Length")
            self.amplicon.grid(row = 1, column = 0,pady=2.5)
            self.amplicon.insert("0", "0")
            self.ampl_label2 = ctk.CTkLabel(self.settings_tabview.tab("Settings"), text= "Base Pairs")
            self.ampl_label2.grid(row = 1, column = 1, padx=10,pady=2.5)
            
            self.DNAconc_label = ctk.CTkLabel(self.settings_tabview.tab("Settings"), text= "DNA Concentration")
            self.DNAconc_label.grid(row = 2, column = 0, pady=2.5)
            self.DNAconc = ctk.CTkEntry(self.settings_tabview.tab("Settings"), placeholder_text="DNA Concentration")
            self.DNAconc.grid(row = 3, column = 0, pady=2.5)
            self.DNAconc.insert("0", "50")
            self.DNAconc_label2 = ctk.CTkLabel(self.settings_tabview.tab("Settings"), text= "mMol")
            self.DNAconc_label2.grid(row = 3, column = 1,padx=10,pady=2.5)
            
            self.temp_label = ctk.CTkLabel(self.settings_tabview.tab("Settings"), text= "Temperature")
            self.temp_label.grid(row = 4, column = 0,pady=2.5)
            self.temp = ctk.CTkEntry(self.settings_tabview.tab("Settings"), placeholder_text="Temperature")
            self.temp.grid(row = 5, column = 0,pady=2.5)
            self.temp.insert("0", "50")
            self.temp_label2 = ctk.CTkLabel(self.settings_tabview.tab("Settings"), text= "C°")
            self.temp_label2.grid(row = 5, column = 1,padx=10,pady=2.5)

            ###ADVANCED SETTINGS###
            self.pc_label = ctk.CTkLabel(self.settings_tabview.tab("Advanced Settings"), text= "Primer Concentration")
            self.pc_label.grid(row = 0, column = 0,pady=2.5)
            self.pc = ctk.CTkEntry(self.settings_tabview.tab("Advanced Settings"), placeholder_text="Primer Concentration")
            self.pc.grid(row = 1, column = 0,pady=2.5)
            self.pc.insert("0", "50")
            self.pc_label2 = ctk.CTkLabel(self.settings_tabview.tab("Advanced Settings"), text= "nMol")
            self.pc_label2.grid(row = 1, column = 1, padx=10,pady=2.5)
            
            self.mc_label = ctk.CTkLabel(self.settings_tabview.tab("Advanced Settings"), text= "Monovalent Ion Concentration")
            self.mc_label.grid(row = 2, column = 0, pady=2.5)
            self.mc = ctk.CTkEntry(self.settings_tabview.tab("Advanced Settings"), placeholder_text="Monovalent Ion Concentration")
            self.mc.grid(row = 3, column = 0, pady=2.5)
            self.mc.insert("0", "50")
            self.mc_label2 = ctk.CTkLabel(self.settings_tabview.tab("Advanced Settings"), text= "mMol")
            self.mc_label2.grid(row = 3, column = 1,padx=10,pady=2.5)

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