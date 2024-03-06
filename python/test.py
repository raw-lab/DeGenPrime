from tkinter import *
from tkinter.filedialog import askopenfile
import customtkinter as ctk
import os, subprocess, webbrowser
from CTkMessagebox import CTkMessagebox

ctk.set_appearance_mode("System")
ctk.set_default_color_theme("dark-blue.json")

class ToplevelWindow(ctk.CTkToplevel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.geometry("400x300")
        self.title("DeGenPrime-Ez Primer Testing Tool")

        def run_test():
            pass

        self.lbl_greeting = ctk.CTkLabel(self,text="Our primer testing tool is compatible with forward and reverse primers. Simply enter your primer sequence and "
                                         + "see if your primer would be filtered by our refinement algorithms.",wraplength = 275, justify= "center")
        self.lbl_greeting.pack(pady=20)
        self.lbl_greeting.pack(padx=30, pady=(25, 5))
        self.select_greeting = ctk.CTkEntry(self, placeholder_text="Enter primer here",justify="center", width = 275)
        self.select_greeting.pack(padx=30, pady=(5, 10))
        self.bp_btn = ctk.CTkButton(self,text="Run Primer Test", command=run_test)
        self.bp_btn.pack(padx=30, pady=10)


class MyTabView(ctk.CTkTabview):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        self.parent = master
        # create tabs
        self.add("Settings")
        self.add("Advanced Settings")
        #self.grid_columnconfigure((0,1,2), weight=1)


        ##old switch for changing measurement, may need later##
        #def switch_event():
        #    pass
            
        #measurements_label = ctk.CTkLabel(master=self.tab("Settings"), text = "Measurement Units")
        #measurements_label.grid(row=0, column=1, padx=(0,5))
        #switch_var = ctk.StringVar(value="on")
        #switch = ctk.CTkSwitch(master=self.tab("Settings"), text="", command=switch_event,
        #                                variable=switch_var, onvalue="on", offvalue="off", fg_color="#939BA2")
        #switch.grid(row=1, column=1)
        #amp_measure_label = ctk.CTkLabel(master=self.tab("Settings"), text = "Amplicon Length")
        #amp_measure_label.grid(row=1, column=0, padx=5)
        #bp_measure_label = ctk.CTkLabel(master=self.tab("Settings"), text = "Base Pairs")
        #bp_measure_label.grid(row=1, column=2, padx=5)
        
        ##SETTINGS##
        min_amp_label = ctk.CTkLabel(master=self.tab("Settings"), text = "Max. Amplicon Length")
        min_amp_label.grid(row=2, column=1, padx=20, pady=10)
        amplicon = ctk.CTkEntry(master=self.tab("Settings"), placeholder_text="Min. Amplicon Length", width = 200, height= 35)
        amplicon.grid(row = 3, column = 1,padx=10,pady=2.5)
        amplicon.insert("0", "100")
        ampl_label2 = ctk.CTkLabel(master=self.tab("Settings"), text= "Base Pairs")
        ampl_label2.grid(row = 3, column = 2, padx=10,pady=2.5)
        
        bpU_label = ctk.CTkLabel(master=self.tab("Settings"), text= "Base Pair Start Point")
        bpU_label.grid(row = 4, column = 0,pady=(30,2.5))
        upper_limit = ctk.CTkEntry(master=self.tab("Settings"), placeholder_text="Start", width = 200, height= 35)
        upper_limit.grid(row = 5, column = 0)
        bpU_label2 = ctk.CTkLabel(master=self.tab("Settings"), text= "to")
        bpU_label2.grid(row = 5, column = 1,pady=(0,10))
        bpL_label = ctk.CTkLabel(master=self.tab("Settings"), text= "Base Pair End Point")
        bpL_label.grid(row = 4, column = 2,pady=(30,2.5))
        lower_limit = ctk.CTkEntry(master=self.tab("Settings"), placeholder_text="End", width = 200, height= 35)
        lower_limit.grid(row = 5, column = 2,pady=2.5)

        templow_label = ctk.CTkLabel(master=self.tab("Settings"), text= "Temperature Lower Bound")
        templow_label.grid(row = 6, column = 0,pady=(30,2.5))
        templow = ctk.CTkEntry(master=self.tab("Settings"), placeholder_text="Lower Bound", width = 200, height= 35)
        templow.grid(row = 7, column = 0)
        templow.insert("0", "50")
        temp_label2 = ctk.CTkLabel(master=self.tab("Settings"), text= "to")
        temp_label2.grid(row = 7, column = 1,pady=(0,10))
        temphigh_label = ctk.CTkLabel(master=self.tab("Settings"), text= "Temperature Upper Bound")
        temphigh_label.grid(row = 6, column = 2,pady=(30,2.5))
        temphigh = ctk.CTkEntry(master=self.tab("Settings"), placeholder_text="Upper Bound", width = 200, height= 35)
        temphigh.grid(row = 7, column = 2,pady=2.5)
        temphigh.insert("0", "65")

        primers_label = ctk.CTkLabel(master=self.tab("Settings"), text= "Preferred # of Primers")
        primers_label.grid(row = 8, column = 1, pady=(30,2.5))
        primers = ctk.CTkEntry(master=self.tab("Settings"), placeholder_text="Number of Primers", width = 200, height= 35)
        primers.grid(row = 9, column = 1, pady=2.5)
        primers.insert("0", "5")
        primers_label2 = ctk.CTkLabel(master=self.tab("Settings"), text= "Primers")
        primers_label2.grid(row = 9, column = 2,padx=10,pady=2.5)

        priMin_label = ctk.CTkLabel(master=self.tab("Settings"), text= "Minimum Primer Length")
        priMin_label.grid(row = 10, column = 0, pady=(30,2.5))
        priMin = ctk.CTkEntry(master=self.tab("Settings"), placeholder_text="Minimum Length", width = 200, height= 35)
        priMin.grid(row = 11, column = 0, pady=2.5)
        priMin.insert("0", "18")
        priMin_label2 = ctk.CTkLabel(master=self.tab("Settings"), text= "to")
        priMin_label2.grid(row = 11, column = 1,pady=(0,20))
        priMax_label = ctk.CTkLabel(master=self.tab("Settings"), text= "Maximum Primer Length")
        priMax_label.grid(row = 10, column = 2, pady=(30,2.5))
        priMax = ctk.CTkEntry(master=self.tab("Settings"), placeholder_text="Maximum Length", width = 200, height= 35)
        priMax.grid(row = 11, column = 2, pady=2.5)
        priMax.insert("0", "25")
        

        ##ADVANCED SETTINGS##
        pc_label = ctk.CTkLabel(master=self.tab("Advanced Settings"), text= "Primer Conc.")
        pc_label.grid(row = 0, column = 1,padx = 20, pady=(10,2.5))
        pc = ctk.CTkEntry(master=self.tab("Advanced Settings"), placeholder_text="Primer Conc.", width = 200, height= 35)
        pc.grid(row = 1, column = 1,padx = 10, pady=10)
        pc.insert("0", "50")
        pc_label2 = ctk.CTkLabel(master=self.tab("Advanced Settings"), text= "nMol")
        pc_label2.grid(row = 1, column = 2, padx=10, pady=10)
        
        mc_label = ctk.CTkLabel(master=self.tab("Advanced Settings"), text= "Monovalent Ion Conc.")
        mc_label.grid(row = 2, column = 1, padx = 20, pady=(10, 2.5))
        mc = ctk.CTkEntry(master=self.tab("Advanced Settings"), placeholder_text="Monovalent Ion Conc.", width = 200, height= 35)
        mc.grid(row = 3, column = 1,padx = 20, pady=10)
        mc.insert("0", "50")
        justify = ctk.CTkLabel(master=self.tab("Advanced Settings"), text= "                                             ")
        justify.grid(row = 3, column = 0, pady=(10,2.5))
        mc_label2 = ctk.CTkLabel(master=self.tab("Advanced Settings"), text= "mMol")
        mc_label2.grid(row = 3, column = 2,padx=10,pady=10)

        #test_lbl = ctk.CTkLabel(master=self.tab("Advanced Settings"), text= " ")
        #test_lbl.grid(row = 4, column = 0, pady=(10,2.5), width = 200, height = 35)
        deltag_label = ctk.CTkLabel(master=self.tab("Advanced Settings"), text= "Gibbs Free Energy")
        deltag_label.grid(row = 4, column = 1, padx = 20, pady=(10, 2.5))
        deltag = ctk.CTkEntry(master=self.tab("Advanced Settings"), placeholder_text="Gibbs Free Energy", width = 200, height= 35)
        deltag.grid(row = 5, column = 1, pady=10)
        deltag.insert("0", "-4")
        deltag_label2 = ctk.CTkLabel(master=self.tab("Advanced Settings"), text= "kCal/mol")
        deltag_label2.grid(row = 5, column = 2,padx = 20, pady=10)



    def pull_settings(self,inputpath):
        def handle_exit_code(exit_code):
                if exit_code == 0:
                    print("Success")
                elif exit_code == 1:
                    CTkMessagebox(title="ERROR", message="ERROR - DeGenPrime has determined your input file is misaligned. Please try a different file.", icon="cancel")
                elif exit_code == 2:
                    CTkMessagebox(title="ERROR", message="ERROR - Settings file not found. Please consider uninstalling and reinstalling DeGenPrime.", icon="cancel")
                elif exit_code == 3:
                    CTkMessagebox(title="ERROR", message="ERROR - Insufficient primers found. Please consider changing settings and re-running.", icon="cancel")
                elif exit_code == 4:
                    CTkMessagebox(title="ERROR", message="ERROR - Invalid input file format. Please try a different file.", icon="cancel")
                elif exit_code == 5:
                    CTkMessagebox(title="ERROR", message="ERROR - Invalid input file. Please try a different file.", icon="cancel")
        #pulling settings from the boxes
        temphigh_get = self.temphigh.get()
        templow_get = self.templow.get()
        pc_get = self.pc.get()
        mc_get = self.mc.get()
        primers_get = self.primers.get()
        amplicon_get = self.amplicon.get()
        priMin_get = self.priMin.get()
        priMax_get = self.priMax.get()
        deltag_get = self.deltag.get()
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
                f"--dG:{deltag_get}",
                inputpath
            ]
        try:
            result = subprocess.run(command, check=True)
            handle_exit_code(result.returncode)
        except subprocess.CalledProcessError as e:
            handle_exit_code(e.returncode)
    


class App(ctk.CTk):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.geometry("950x600")
        self.title("DeGenPrime-Ez")

    
        def find_input_file():
            file = askopenfile(mode ='r', filetypes =[('Acceptable Filetypes', '*.fasta *.fna *.clust *.faa')])
            if file:
                inputpath = os.path.abspath(file.name)
                self.lbl_inputpath = ctk.CTkLabel(self, text="Your selected file is: " + str(inputpath),wraplength = 300,justify="center")
                self.lbl_inputpath.grid(row=2, column=2, padx = 10, pady=10)
        def open_browser():
            url = 'https://github.com/raw-lab/DeGenPrime'
            webbrowser.open_new(url)

        self.sidebar = ctk.CTkFrame(self, width=250)
        self.sidebar.pack(side="left", fill="y")

        self.tab_view = MyTabView(self)
        self.tab_view.pack(side = "right", fill = "both", expand = True, padx=20, pady=(10,20))

        self.test_primers = ctk.CTkButton(self.sidebar, text="Test Primers", command=self.open_toplevel, height = 50)
        self.test_primers.pack(padx=20, pady=50)

        self.browse_files = ctk.CTkButton(self.sidebar, text="Browse Files", command=find_input_file, height = 50)
        self.browse_files.pack(padx=20, pady=50)

        self.run_program = ctk.CTkButton(self.sidebar, text="Run Program", command= self.tab_view.pull_settings, height = 50)
        self.run_program.pack(padx=20, pady=50)

        self.github = ctk.CTkButton(self.sidebar, text="GitHub", command=open_browser, height = 50)
        self.github.pack(padx=20, pady=50)

        self.toplevel_window = None

    def open_toplevel(self):
        if self.toplevel_window is None or not self.toplevel_window.winfo_exists():
            self.toplevel_window = ToplevelWindow(self)
        else:
            self.toplevel_window.focus()


app = App()
app.mainloop()