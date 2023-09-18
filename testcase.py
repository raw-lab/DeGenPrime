import customtkinter as ctk

def on_checkbox_click(checkbox):
    global checkbox1_state, checkbox2_state

    if checkbox == 1:
        if (checkbox1_state == 0):
            checkbox2.configure(state = 'disabled', text = checkbox1_state) #checkbox 2 is disabled, text changes
        else:
            checkbox2.configure(state = 'enabled', text = checkbox1_state)
    elif checkbox == 2:
        checkbox2_state = not checkbox2_state
        checkbox1.configure(state= 'disabled', text = checkbox2_state)

root = ctk.CTk()

checkbox1_state = False
checkbox2_state = False

checkbox1 = ctk.CTkCheckBox(root, text="Option 1", command=lambda: on_checkbox_click(1))
checkbox2 = ctk.CTkCheckBox(root, text="Option 2", command=lambda: on_checkbox_click(2))

checkbox1.pack(pady=10)
checkbox2.pack(pady=10)

root.mainloop()