from tkinter import *
from tkinter.filedialog import askopenfile
import customtkinter as ctk
import os, subprocess

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
        def amp_click():
            
            # configure grid layout (4x4)
            self.grid_columnconfigure((0,1,2,3), weight=1)
            self.grid_rowconfigure((0, 1, 2), weight=1)

            #delete old info
            self.bp_btn.destroy()
            self.amp_btn.destroy()
            self.lbl_greeting.destroy()
            self.select_greeting.destroy()
            
            #event handling for finding the input file
            def find_input_file():
                file = askopenfile(mode ='r', filetypes =[('Acceptable Filetypes', '*.fasta *.fna *.clust *.faa')])
                if file:
                    global inputpath
                    inputpath = os.path.abspath(file.name)
                    self.lbl_inputpath = ctk.CTkLabel(self, text="Your selected file is: " + str(inputpath),wraplength = 300,justify="center")
                    self.lbl_inputpath.grid(row=2, column=2, padx = 10, pady=10)
                    filecontent = file.read()

            self.grid_rowconfigure((0,1), weight=1)
            self.grid_columnconfigure((0, 1, 2), weight=1)

            #initializing the scrollbar
            self.settings_scroll = ctk.CTkScrollableFrame(self)
            self.settings_scroll.grid(row=1, column=2, sticky="nsew")
            self.settings_scroll.grid_columnconfigure(0, weight=1)
            
            self.ampl_label = ctk.CTkLabel(self.settings_scroll, text= "Min. Amplicon Length")
            self.ampl_label.grid(row = 2, column = 0,pady=2.5)
            self.amplicon = ctk.CTkEntry(self.settings_scroll, placeholder_text="Min. Amplicon Length")
            self.amplicon.grid(row = 3, column = 0,pady=2.5)
            self.amplicon.insert("0", "1")
            self.ampl_label2 = ctk.CTkLabel(self.settings_scroll, text= "Base Pairs")
            self.ampl_label2.grid(row = 3, column = 1, padx=10,pady=2.5)
            
            #create the frame for the three side buttons -- browse, run, and save
            self.sidebar_frame = ctk.CTkFrame(self, width=140, corner_radius=0)
            self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
            self.sidebar_frame.grid_rowconfigure(4, weight=1)
            self.btn_browse = ctk.CTkButton(self.sidebar_frame, text = "BROWSE FILES", command = find_input_file)
            self.btn_browse.grid(row=0, column=0, padx=20, pady=(40,30))
            self.btn_run = ctk.CTkButton(self.sidebar_frame, text = "RUN PROGRAM", command = run_amp)
            self.btn_run.grid(row=2, column=0, padx=20, pady=30)
            self.btn_save = ctk.CTkButton(self.sidebar_frame, text = "PLACEHOLDER")
            self.btn_save.grid(row=3, column=0, padx=20, pady=30)


            #upper bound of temperature
            self.temphigh_label = ctk.CTkLabel(self.settings_scroll, text= "Temp. (Upper Bound)",justify="center")
            self.temphigh_label.grid(row = 6, column = 0, pady=2.5)
            self.temphigh = ctk.CTkEntry(self.settings_scroll, placeholder_text="DNA Temperature (Upper Bound)")
            self.temphigh.grid(row = 7, column = 0, pady=2.5)
            self.temphigh.insert("0", "65")
            self.temphigh_label2 = ctk.CTkLabel(self.settings_scroll, text= "C°")
            self.temphigh_label2.grid(row = 7, column = 1,padx=10,pady=2.5)
            
            #lower bound of temperature
            self.templow_label = ctk.CTkLabel(self.settings_scroll, text= "Temp. (Lower Bound)")
            self.templow_label.grid(row = 8, column = 0,pady=2.5)
            self.templow = ctk.CTkEntry(self.settings_scroll, placeholder_text="Temperature (Lower Bound)")
            self.templow.grid(row = 9, column = 0,pady=2.5)
            self.templow.insert("0", "50")
            self.templow_label2 = ctk.CTkLabel(self.settings_scroll, text= "C°")
            self.templow_label2.grid(row = 9, column = 1,padx=10,pady=2.5)

            #primer concentration
            self.pc_label = ctk.CTkLabel(self.settings_scroll, text= "Primer Conc.")
            self.pc_label.grid(row = 10, column = 0,pady=2.5)
            self.pc = ctk.CTkEntry(self.settings_scroll, placeholder_text="Primer Conc.")
            self.pc.grid(row = 11, column = 0,pady=2.5)
            self.pc.insert("0", "50")
            self.pc_label2 = ctk.CTkLabel(self.settings_scroll, text= "nMol")
            self.pc_label2.grid(row = 11, column = 1, padx=10,pady=2.5)
            
            #monovalent ion/salt concentration
            self.mc_label = ctk.CTkLabel(self.settings_scroll, text= "Monovalent Ion Conc.")
            self.mc_label.grid(row = 12, column = 0, pady=2.5)
            self.mc = ctk.CTkEntry(self.settings_scroll, placeholder_text="Monovalent Ion Conc.")
            self.mc.grid(row = 13, column = 0, pady=2.5)
            self.mc.insert("0", "50")
            self.mc_label2 = ctk.CTkLabel(self.settings_scroll, text= "mMol")
            self.mc_label2.grid(row = 13, column = 1,padx=10,pady=2.5)

            #preferred amount of primers
            self.primers_label = ctk.CTkLabel(self.settings_scroll, text= "Preferred # of Primers")
            self.primers_label.grid(row = 14, column = 0, pady=2.5)
            self.primers = ctk.CTkEntry(self.settings_scroll, placeholder_text="Number of Primers")
            self.primers.grid(row = 15, column = 0, pady=2.5)
            self.primers.insert("0", "5")
            self.primers_label2 = ctk.CTkLabel(self.settings_scroll, text= "Primers")
            self.primers_label2.grid(row = 15, column = 1,padx=10,pady=2.5)

            #preferred amount of primers
            self.priMin_label = ctk.CTkLabel(self.settings_scroll, text= "Minimum Length of Primers")
            self.priMin_label.grid(row = 16, column = 0, pady=2.5)
            self.priMin = ctk.CTkEntry(self.settings_scroll, placeholder_text="Minimum Length")
            self.priMin.grid(row = 17, column = 0, pady=2.5)
            self.priMin.insert("0", "18")
            self.priMin_label2 = ctk.CTkLabel(self.settings_scroll, text= "Primers")
            self.priMin_label2.grid(row = 17, column = 1,padx=10,pady=2.5)

            self.priMax_label = ctk.CTkLabel(self.settings_scroll, text= "Minimum Length of Primers")
            self.priMax_label.grid(row = 18, column = 0, pady=2.5)
            self.priMax = ctk.CTkEntry(self.settings_scroll, placeholder_text="Minimum Length")
            self.priMax.grid(row = 19, column = 0, pady=2.5)
            self.priMax.insert("0", "25")
            self.priMax_label2 = ctk.CTkLabel(self.settings_scroll, text= "Primers")
            self.priMax_label2.grid(row = 19, column = 1,padx=10,pady=2.5)
            



        def run_amp():
            #pulling settings from the boxes
            temphigh_get = self.temphigh.get()
            if(int(temphigh_get) > 50): temphigh_get = 50
            templow_get = self.templow.get()
            if(int(templow_get) < 0): templow_get = 50
            pc_get = self.pc.get()
            if(int(pc_get) <= 0): pc_get = 50
            mc_get = self.mc.get()
            if(int(mc_get) <= 0): mc_get = 50
            primers_get = self.primers.get()
            if(int(primers_get) > 10 or int(primers_get) <= 0): primers_get = 10
            amplicon_get = self.amplicon.get()
            if(int(amplicon_get) < 0): amplicon_get = 0
            priMin_get = self.priMin.get()
            if(int(priMin_get) < 18): priMin_get = 20
            priMax_get = self.priMax.get()
            if(int(priMax_get) < 18): priMax_get = 22
            command = [
                    "degenprime",
                    f"--amplicon:{amplicon_get}",
                    f"--min_temp:{templow_get}",
                    f"--max_temp:{temphigh_get}",
                    f"--primer_conc:{pc_get}",
                    f"--salt_conc:{mc_get}",
                    f"--max_primers:{primers_get}",
                    f"--min_primer_len:{priMin_get}",
                    f"--max_primer_len:{priMax_get}",
                    inputpath
                ]
                    #run command, save primer to .txt file
            output = subprocess.check_output(command, text=True)
            output_file_address = "primer.txt"
            with open(output_file_address, 'w') as file:
                file.write(output)
        
        def bp_click():
            
            # configure grid layout (4x4)
            self.grid_columnconfigure((0,1,2,3), weight=1)
            self.grid_rowconfigure((0, 1, 2), weight=1)

            #delete old info
            self.bp_btn.destroy()
            self.amp_btn.destroy()
            self.lbl_greeting.destroy()
            self.select_greeting.destroy()
            
            #event handling for finding the input file
            def find_input_file():
                file = askopenfile(mode ='r', filetypes =[('Acceptable Filetypes', '*.fasta *.fna *.clust *.faa')])
                if file:
                    global inputpath
                    inputpath = os.path.abspath(file.name)
                    self.lbl_inputpath = ctk.CTkLabel(self, text="Your selected file is: " + str(inputpath),wraplength = 300,justify="center")
                    self.lbl_inputpath.grid(row=2, column=2, padx = 10, pady=10)
                    filecontent = file.read()

            self.grid_rowconfigure((0,1), weight=1)
            self.grid_columnconfigure((0, 1, 2), weight=1)

            #initializing the scrollbar
            self.settings_scroll = ctk.CTkScrollableFrame(self)
            self.settings_scroll.grid(row=1, column=2, sticky="nsew")
            self.settings_scroll.grid_columnconfigure(0, weight=1)
            
            self.bpU_label = ctk.CTkLabel(self.settings_scroll, text= "Base Pair Upper Limit")
            self.bpU_label.grid(row = 2, column = 0,pady=2.5)
            self.upper_limit = ctk.CTkEntry(self.settings_scroll, placeholder_text="Upper Limit")
            self.upper_limit.grid(row = 3, column = 0,pady=2.5)
            self.bpU_label2 = ctk.CTkLabel(self.settings_scroll, text= "Position")
            self.bpU_label2.grid(row = 3, column = 1, padx=10,pady=2.5)
            self.bpL_label = ctk.CTkLabel(self.settings_scroll, text= "Base Pair Lower Limit")
            self.bpL_label.grid(row = 4, column = 0,pady=2.5)
            self.lower_limit = ctk.CTkEntry(self.settings_scroll, placeholder_text="Lower Limit")
            self.lower_limit.grid(row = 5, column = 0,pady=2.5)
            self.bpL_label2 = ctk.CTkLabel(self.settings_scroll, text= "Position")
            self.bpL_label2.grid(row = 5, column = 1, padx=10,pady=2.5)
            
            #create the frame for the three side buttons -- browse, run, and save
            self.sidebar_frame = ctk.CTkFrame(self, width=140, corner_radius=0)
            self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
            self.sidebar_frame.grid_rowconfigure(4, weight=1)
            self.btn_browse = ctk.CTkButton(self.sidebar_frame, text = "BROWSE FILES", command = find_input_file)
            self.btn_browse.grid(row=0, column=0, padx=20, pady=(40,30))
            self.btn_run = ctk.CTkButton(self.sidebar_frame, text = "RUN PROGRAM", command = run_basepair)
            self.btn_run.grid(row=2, column=0, padx=20, pady=30)
            self.btn_save = ctk.CTkButton(self.sidebar_frame, text = "PLACEHOLDER")
            self.btn_save.grid(row=3, column=0, padx=20, pady=30)


            #upper bound of temperature
            self.temphigh_label = ctk.CTkLabel(self.settings_scroll, text= "Temp. (Upper Bound)",justify="center")
            self.temphigh_label.grid(row = 6, column = 0, pady=2.5)
            self.temphigh = ctk.CTkEntry(self.settings_scroll, placeholder_text="DNA Temperature (Upper Bound)")
            self.temphigh.grid(row = 7, column = 0, pady=2.5)
            self.temphigh.insert("0", "65")
            self.temphigh_label2 = ctk.CTkLabel(self.settings_scroll, text= "C°")
            self.temphigh_label2.grid(row = 7, column = 1,padx=10,pady=2.5)
            
            #lower bound of temperature
            self.templow_label = ctk.CTkLabel(self.settings_scroll, text= "Temp. (Lower Bound)")
            self.templow_label.grid(row = 8, column = 0,pady=2.5)
            self.templow = ctk.CTkEntry(self.settings_scroll, placeholder_text="Temperature (Lower Bound)")
            self.templow.grid(row = 9, column = 0,pady=2.5)
            self.templow.insert("0", "50")
            self.templow_label2 = ctk.CTkLabel(self.settings_scroll, text= "C°")
            self.templow_label2.grid(row = 9, column = 1,padx=10,pady=2.5)

            #primer concentration
            self.pc_label = ctk.CTkLabel(self.settings_scroll, text= "Primer Conc.")
            self.pc_label.grid(row = 10, column = 0,pady=2.5)
            self.pc = ctk.CTkEntry(self.settings_scroll, placeholder_text="Primer Conc.")
            self.pc.grid(row = 11, column = 0,pady=2.5)
            self.pc.insert("0", "50")
            self.pc_label2 = ctk.CTkLabel(self.settings_scroll, text= "nMol")
            self.pc_label2.grid(row = 11, column = 1, padx=10,pady=2.5)
            
            #monovalent ion/salt concentration
            self.mc_label = ctk.CTkLabel(self.settings_scroll, text= "Monovalent Ion Conc.")
            self.mc_label.grid(row = 12, column = 0, pady=2.5)
            self.mc = ctk.CTkEntry(self.settings_scroll, placeholder_text="Monovalent Ion Conc.")
            self.mc.grid(row = 13, column = 0, pady=2.5)
            self.mc.insert("0", "50")
            self.mc_label2 = ctk.CTkLabel(self.settings_scroll, text= "mMol")
            self.mc_label2.grid(row = 13, column = 1,padx=10,pady=2.5)

            #preferred amount of primers
            self.primers_label = ctk.CTkLabel(self.settings_scroll, text= "Preferred # of Primers")
            self.primers_label.grid(row = 14, column = 0, pady=2.5)
            self.primers = ctk.CTkEntry(self.settings_scroll, placeholder_text="Number of Primers")
            self.primers.grid(row = 15, column = 0, pady=2.5)
            self.primers.insert("0", "5")
            self.primers_label2 = ctk.CTkLabel(self.settings_scroll, text= "Primers")
            self.primers_label2.grid(row = 15, column = 1,padx=10,pady=2.5)

            #preferred length of primers
            self.priMin_label = ctk.CTkLabel(self.settings_scroll, text= "Minimum Length of Primers")
            self.priMin_label.grid(row = 16, column = 0, pady=2.5)
            self.priMin = ctk.CTkEntry(self.settings_scroll, placeholder_text="Minimum Length")
            self.priMin.grid(row = 17, column = 0, pady=2.5)
            self.priMin.insert("0", "18")
            self.priMin_label2 = ctk.CTkLabel(self.settings_scroll, text= "Primers")
            self.priMin_label2.grid(row = 17, column = 1,padx=10,pady=2.5)

            self.priMax_label = ctk.CTkLabel(self.settings_scroll, text= "Minimum Length of Primers")
            self.priMax_label.grid(row = 18, column = 0, pady=2.5)
            self.priMax = ctk.CTkEntry(self.settings_scroll, placeholder_text="Minimum Length")
            self.priMax.grid(row = 19, column = 0, pady=2.5)
            self.priMax.insert("0", "25")
            self.priMax_label2 = ctk.CTkLabel(self.settings_scroll, text= "Primers")
            self.priMax_label2.grid(row = 19, column = 1,padx=10,pady=2.5)


        def run_basepair():
            #pulling settings from the boxes
            temphigh_get = self.temphigh.get()
            if(int(temphigh_get) > 50): temphigh_get = 50
            templow_get = self.templow.get()
            if(int(templow_get) < 0): templow_get = 50
            pc_get = self.pc.get()
            if(int(pc_get) <= 0): pc_get = 50
            mc_get = self.mc.get()
            if(int(mc_get) <= 0): mc_get = 50
            primers_get = self.primers.get()
            if(int(primers_get) > 10 or int(primers_get) <= 0): primers_get = 10
            bpU_get = self.upper_limit.get()
            bpL_get = self.lower_limit.get()
            if(int(bpL_get) < 0): bpL_get = 0
            priMin_get = self.priMin.get()
            if(int(priMin_get) < 18): priMin_get = 20
            priMax_get = self.priMax.get()
            if(int(priMax_get) < 18): priMax_get = 22
            command = [
                    "degenprime",
                    f"--begin:{bpL_get}",
                    f"--end:{bpU_get}",
                    f"--min_temp:{templow_get}",
                    f"--max_temp:{temphigh_get}",
                    f"--primer_conc:{pc_get}",
                    f"--salt_conc:{mc_get}",
                    f"--max_primers:{primers_get}",
                    f"--min_primer_len:{priMin_get}",
                    f"--max_primer_len:{priMax_get}",
                    inputpath
                ]
            output = subprocess.check_output(command, text=True)
            output_file_address = "primer.txt"
            with open(output_file_address, 'w') as file:
                file.write(output)


        #original welcome message frame
        self.sidebar_frame = ctk.CTkFrame(self, width = 500)
        self.lbl_greeting = ctk.CTkLabel(self,text="Welcome to DeGenPrime. Our goal is to maximize your primer design experience while minimizing your headaches. "
                                         + "Input either an unaligned or a prealigned .fasta or .fna file, and allow DeGenPrime to recommend the best primers for your PCR use."
                                          ,wraplength = 450, justify= "center")
        self.lbl_greeting.grid(padx=30, pady=(50, 5))
        self.select_greeting = ctk.CTkLabel(self, text="Please select how you would like to measure your primer. To change your selection, restart the program.",wraplength=450,justify="center")
        self.select_greeting.grid(padx=30, pady=(5, 10))
        self.amp_btn = ctk.CTkButton(self,text="Amplicon Length", command = amp_click)
        self.amp_btn.grid(padx=30, pady=10)
        self.bp_btn = ctk.CTkButton(self,text="Base Pair Ranges", command= bp_click)
        self.bp_btn.grid(padx=30, pady=10)

#initializes program
if __name__ == "__main__":
    app = App()
    app.mainloop()