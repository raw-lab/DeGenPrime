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
            
            #creates the toggle switch function
            def toggle_switch(choice):
                if(choice == "Amplicon Length"):
                    self.ampl_label = ctk.CTkLabel(self.settings_scroll, text= "Min. Amplicon Length")
                    self.ampl_label.grid(row = 2, column = 0,pady=2.5)
                    self.amplicon = ctk.CTkEntry(self.settings_scroll, placeholder_text="Min. Amplicon Length")
                    self.amplicon.grid(row = 3, column = 0,pady=2.5)
                    self.amplicon.insert("0", "1")
                    self.ampl_label2 = ctk.CTkLabel(self.settings_scroll, text= "Base Pairs")
                    self.ampl_label2.grid(row = 3, column = 1, padx=10,pady=2.5)
                    self.bpU_label.destroy()
                    self.upper_limit.destroy()
                    self.bpU_label2.destroy()
                    self.bpL_label.destroy()
                    self.lower_limit.destroy()
                    self.bpL_label2.destroy()
                elif(choice == "Base Pairs"):
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
                    self.ampl_label.destroy()
                    self.amplicon.destroy()
                    self.ampl_label2.destroy()

            #toggle switch from amplicon length and range of base pairs
            self.switcher_label = ctk.CTkLabel(self.settings_scroll, text="Measurement Type")
            self.switcher_label.grid(row=0, column=0,padx=3, pady= 2.5)
            self.length_switcher = ctk.CTkComboBox(self.settings_scroll, values=["Amplicon Length","Base Pairs"], command=toggle_switch)
            self.length_switcher.grid(row=1, column=0,padx=10, pady=2.5)

            #legacy code for saving settings to a .csv - DO NOT USE.
            #def save_settings():
                #getting our amplicon length or base pair range based on global choice state
                #choice_state = str(self.length_switcher.get())
                #if(choice_state == "Amplicon Length"):
                #    amplicon_get = self.amplicon.get()
                #    if(int(amplicon_get) < 0): amplicon_get = 0
                #elif(choice_state == "Base Pairs"):
                #    bpU_get = self.upper_limit.get()
                #    bpL_get = self.lower_limit.get()
                #    if(int(bpL_get) < 0): bpL_get = 0
                #temphigh_get = self.temphigh.get()
                #if(int(temphigh_get) < 0): temphigh_get = 50
                #templow_get = self.templow.get()
                #if(int(templow_get) < 50): templow_get = 50
                #pc_get = self.pc.get()
                #if(int(pc_get) <= 0): pc_get = 50
                #mc_get = self.mc.get()
                #if(int(mc_get) <= 0): mc_get = 50

                #with open('SETTINGS.csv', 'w', newline='') as file:
                #    writer = csv.writer(file)
                #    if(choice_state == "Amplicon Length"):
                #        field = ["Amplicon_Length", "Temperature_High", "Temperature_Low", "Primer_Concentration", "Monovalent_Ion_Concentration"]
                #        writer.writerow(field)
                #        writer.writerow([amplicon_get, temphigh_get, templow_get, pc_get, mc_get])
                #    elif(choice_state == "Base Pairs"):
                #        field = ["BasePair_Upper", "BasePair_Lower", "Temperature_High", "Temperature_Low", "Primer_Concentration", "Monovalent_Ion_Concentration"]
                #        writer.writerow(field)
                #        writer.writerow([bpU_get, bpL_get, temphigh_get, templow_get, pc_get, mc_get])
            
            #create the frame for the three side buttons -- browse, run, and save
            self.sidebar_frame = ctk.CTkFrame(self, width=140, corner_radius=0)
            self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
            self.sidebar_frame.grid_rowconfigure(4, weight=1)
            self.btn_browse = ctk.CTkButton(self.sidebar_frame, text = "BROWSE FILES", command = find_input_file)
            self.btn_browse.grid(row=0, column=0, padx=20, pady=(40,30))
            self.btn_run = ctk.CTkButton(self.sidebar_frame, text = "RUN PROGRAM", command = run_program)
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



        def run_program():
            #pulling settings from the boxes
            temphigh_get = self.temphigh.get()
            if(int(temphigh_get) < 0): temphigh_get = 50
            templow_get = self.templow.get()
            if(int(templow_get) < 50): templow_get = 50
            pc_get = self.pc.get()
            if(int(pc_get) <= 0): pc_get = 50
            mc_get = self.mc.get()
            if(int(mc_get) <= 0): mc_get = 50
            choice_state = str(self.length_switcher.get())
            primers_get = self.primers.get()
            if(int(primers_get) > 10 or int(primers_get) <= 0): primers_get = 10
            if(choice_state == "Amplicon Length"):
                    amplicon_get = self.amplicon.get()
                    if(int(amplicon_get) < 0): amplicon_get = 0
                    command = [
                         "./degenprime",
                         f"--amplicon:{amplicon_get}",
                         f"--min_temp:{templow_get}",
                         f"--max_temp:{temphigh_get}",
                         f"--primer_conc:{pc_get}",
                         f"--salt_conc:{mc_get}",
                         f"--max_primers:{primers_get}",
                         str(inputpath)
                        ]
                    #run command, save primer to .txt file
                    output = subprocess.check_output(command, text=True)
                    output_file_address = "primer.txt"
                    with open(output_file_address, 'w') as file:
                        file.write(output)
            elif(choice_state == "Base Pairs"):
                    bpU_get = self.upper_limit.get()
                    bpL_get = self.lower_limit.get()
                    if(int(bpL_get) < 0): bpL_get = 0
                    command = [
                         "./degenprime",
                         f"--begin:{bpL_get}",
                         f"--end:{bpU_get}",
                         f"--min_temp:{templow_get}",
                         f"--max_temp:{temphigh_get}",
                         f"--primer_conc:{pc_get}",
                         f"--salt_conc:{mc_get}",
                         f"--max_primers:{primers_get}",
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
        self.lbl_greeting.grid(padx=30, pady=50)
        self.btn_greeting = ctk.CTkButton(self,text="ENTER")
        self.btn_greeting.bind("<Button-1>", handle_click)
        self.btn_greeting.grid(padx=30, pady=10)

#initializes program
if __name__ == "__main__":
    app = App()
    app.mainloop()