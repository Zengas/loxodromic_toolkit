import tkinter as tk
from tkinter import ttk, messagebox
import math
from rhumb import Rhumb

def build_gui(parent):
    def ddm_to_decimal(degrees, minutes, hemi):
        degrees = float(degrees)
        minutes = float(minutes)
        value = degrees + minutes / 60
        if hemi in ("S", "W"):
            value = -value
        return value

    def decimal_to_ddm(val, hemi_pos, hemi_neg, deg_digits=2):
        hemi = hemi_pos if val >= 0 else hemi_neg
        abs_val = abs(val)
        degrees = int(abs_val)
        minutes = (abs_val - degrees) * 60
        deg_str = str(degrees).zfill(deg_digits)
        min_str = f"{minutes:05.2f}"
        return f"{deg_str}° {min_str}' {hemi}"

    # Storage for latest calculation results (for the graph)
    latest_azimuths = {"alpha1": None, "alpha2": None, "distance_nm": None}

    # ========== Title ==========
    title = ttk.Label(
        parent, 
        text="Enter Departing & Arriving Coordinates (DDM), Choose Segment Count",
        font=("Arial", 10, "italic")
    )
    title.pack(pady=(10, 2))

    main_frame = ttk.Frame(parent)
    main_frame.pack(pady=5, padx=5)
    
    # ========== EXPLANATION ==========
    EXPLANATION_MESSAGE = (
    "\nCalculates the distance between two points of the Earth specified geodesic (geographical) coordinates along the shortest path - the great circle (orthodrome). Calculates the initial and final course angles and azimuth at intermediate points between the two given.."
    )

    # --- Inline Input Frame with Vertical Separator ---
    input_frame = ttk.Frame(main_frame)
    input_frame.pack(pady=(0,8), anchor="center")

    dep_lat_hemi = tk.StringVar(value="N")
    dep_lon_hemi = tk.StringVar(value="W")
    arr_lat_hemi = tk.StringVar(value="N")
    arr_lon_hemi = tk.StringVar(value="W")

    def coordinate_section(parent, label, lat_hemi_var, lon_hemi_var):
        section = ttk.LabelFrame(parent, text=label, padding=10)
        ttk.Label(section, text="Latitude:").grid(row=0, column=0, sticky="e", pady=2)
        lat_deg = ttk.Entry(section, width=6); lat_deg.grid(row=0, column=1)
        ttk.Label(section, text="°").grid(row=0, column=2)
        lat_min = ttk.Entry(section, width=8); lat_min.grid(row=0, column=3)
        ttk.Label(section, text="'").grid(row=0, column=4)
        ttk.Radiobutton(section, text="N", variable=lat_hemi_var, value="N").grid(row=0, column=5)
        ttk.Radiobutton(section, text="S", variable=lat_hemi_var, value="S").grid(row=0, column=6)
        ttk.Label(section, text="Longitude:").grid(row=1, column=0, sticky="e", pady=2)
        lon_deg = ttk.Entry(section, width=6); lon_deg.grid(row=1, column=1)
        ttk.Label(section, text="°").grid(row=1, column=2)
        lon_min = ttk.Entry(section, width=8); lon_min.grid(row=1, column=3)
        ttk.Label(section, text="'").grid(row=1, column=4)
        ttk.Radiobutton(section, text="W", variable=lon_hemi_var, value="W").grid(row=1, column=5)
        ttk.Radiobutton(section, text="E", variable=lon_hemi_var, value="E").grid(row=1, column=6)
        return section, lat_deg, lat_min, lon_deg, lon_min

    dep_section, dep_lat_deg, dep_lat_min, dep_lon_deg, dep_lon_min = coordinate_section(
        input_frame, "Departing Coordinates", dep_lat_hemi, dep_lon_hemi
    )
    arr_section, arr_lat_deg, arr_lat_min, arr_lon_deg, arr_lon_min = coordinate_section(
        input_frame, "Arriving Coordinates", arr_lat_hemi, arr_lon_hemi
    )
    dep_section.pack(side="left", padx=(0, 14))
    ttk.Separator(input_frame, orient="vertical").pack(side="left", fill="y", pady=2)
    arr_section.pack(side="left", padx=(14, 0))

    seg_frame = ttk.Frame(main_frame)
    seg_frame.pack(pady=(0, 8), anchor="center")
    ttk.Label(seg_frame, text="Number of Segments:").pack(side="left", padx=(0, 8))
    segments_var = tk.StringVar(value="10")
    segments_choices = [str(x) for x in range(10, 101, 10)]
    segments_dropdown = ttk.Combobox(seg_frame, textvariable=segments_var, values=segments_choices, width=4, state="readonly")
    segments_dropdown.pack(side="left")

    btn_frame = ttk.Frame(parent)
    btn_frame.pack(pady=(5, 0))

    # --- Result Frame with Scrollbars ---
    result_frame = ttk.Frame(parent)
    result_frame.pack(pady=10, fill="both", expand=True)

    result_text = tk.Text(
        result_frame,
        height=12,
        width=66,
        font=("Courier New", 11),
        wrap="word",
        state="normal"  # Start as normal to insert the message
    )
    result_text.insert(tk.END, EXPLANATION_MESSAGE)
    result_text.config(state="disabled")
    vscroll = ttk.Scrollbar(result_frame, orient="vertical", command=result_text.yview)
    result_text.configure(yscrollcommand=vscroll.set)
    result_text.grid(row=0, column=0, sticky="nsew")
    vscroll.grid(row=0, column=1, sticky="ns")
    result_frame.rowconfigure(0, weight=1)
    result_frame.columnconfigure(0, weight=1)

    def calculate():
        try:
            r = Rhumb()
            lat1 = ddm_to_decimal(dep_lat_deg.get(), dep_lat_min.get(), dep_lat_hemi.get())
            lon1 = ddm_to_decimal(dep_lon_deg.get(), dep_lon_min.get(), dep_lon_hemi.get())
            lat2 = ddm_to_decimal(arr_lat_deg.get(), arr_lat_min.get(), arr_lat_hemi.get())
            lon2 = ddm_to_decimal(arr_lon_deg.get(), arr_lon_min.get(), arr_lon_hemi.get())
            segs = int(segments_var.get())
            s, alpha1, alpha2 = r.vincenty_inverse(lat1, lon1, lat2, lon2)
            distance_nm = s / 1852.0
            segment_len = s / segs
            segment_nm = segment_len / 1852.0
            waypoints = []
            for i in range(segs + 1):
                d = segment_len * i
                phi, lam, az = r.vincenty_direct(lat1, lon1, alpha1, d)
                lat_ddm = decimal_to_ddm(phi, "N", "S", deg_digits=2)
                lon_ddm = decimal_to_ddm(lam, "E", "W", deg_digits=3)
                az_str = f"{az:6.2f}°"
                waypoints.append(f"{i+1:02d}: {az_str}   {lat_ddm}   {lon_ddm}")
            result_str = (
                f"--- Great Circle Calculation Result (WGS84 Orthodrome) ---\n\n"
                f"Initial Azimuth : {alpha1:6.2f}°\n"
                f"Final Azimuth   : {alpha2:6.2f}°\n"
                f"Segment Distance: {segment_nm:,.2f} NM\n"
                f"Total Distance  : {distance_nm:,.2f} NM\n\n"
                f"---------------- Waypoints ----------------\n\n"
                f"    Azimuth   Latitude       Longitude\n"
            )
            result_str += "\n".join(waypoints)
            result_text.config(state="normal")
            result_text.delete("1.0", tk.END)
            result_text.insert(tk.END, result_str)
            result_text.config(state="disabled")
            # Store for graph
            latest_azimuths["alpha1"] = alpha1
            latest_azimuths["alpha2"] = alpha2
            latest_azimuths["distance_nm"] = distance_nm
        except Exception as e:
            messagebox.showerror("Error", str(e), parent=parent)

    def clear():
        for e in [dep_lat_deg, dep_lat_min, dep_lon_deg, dep_lon_min,
                  arr_lat_deg, arr_lat_min, arr_lon_deg, arr_lon_min]:
            e.delete(0, tk.END)
        dep_lat_hemi.set("N")
        dep_lon_hemi.set("W")
        arr_lat_hemi.set("N")
        arr_lon_hemi.set("W")
        segments_var.set("10")
        result_text.config(state="normal")
        result_text.delete("1.0", tk.END)
        result_text.insert(tk.END, EXPLANATION_MESSAGE)   # <-- Insert message again
        result_text.config(state="disabled")
        latest_azimuths["alpha1"] = None
        latest_azimuths["alpha2"] = None
        latest_azimuths["distance_nm"] = None

    def show_graph():
        try:
            alpha1 = latest_azimuths["alpha1"]
            alpha2 = latest_azimuths["alpha2"]
            distance_nm = latest_azimuths["distance_nm"]

            if alpha1 is None or alpha2 is None:
                messagebox.showinfo("Graph", "Please perform a calculation first.", parent=parent)
                return

            popup = tk.Toplevel(parent)
            popup.title("Great Circle Azimuths Graph")
            try:
                popup.iconbitmap("compass256.ico")
            except Exception:
                pass
            popup.geometry("350x410")
            popup.resizable(False, False)

            canvas = tk.Canvas(popup, width=330, height=330, bg="white")
            canvas.pack(pady=5)

            center_x = 165
            center_y = 165
            radius = 120

            # Circle
            canvas.create_oval(center_x - radius, center_y - radius, center_x + radius, center_y + radius, outline="black", width=2)
            font_c = ("Arial", 12, "bold")
            canvas.create_text(center_x, center_y - radius - 15, text="N", font=font_c)
            canvas.create_text(center_x + radius + 15, center_y, text="E", font=font_c)
            canvas.create_text(center_x, center_y + radius + 15, text="S", font=font_c)
            canvas.create_text(center_x - radius - 15, center_y, text="W", font=font_c)
            canvas.create_line(center_x, center_y - radius, center_x, center_y + radius, dash=(2, 2))
            canvas.create_line(center_x - radius, center_y, center_x + radius, center_y, dash=(2, 2))

            arrow_len = radius * 0.92

            # Departing azimuth (blue arrow)
            angle_rad_dep = math.radians(alpha1 - 90)
            x_end_dep = center_x + arrow_len * math.cos(angle_rad_dep)
            y_end_dep = center_y + arrow_len * math.sin(angle_rad_dep)
            canvas.create_line(center_x, center_y, x_end_dep, y_end_dep, fill="blue", width=4, arrow=tk.LAST)

            # Arrival azimuth (red arrow)
            angle_rad_arr = math.radians(alpha2 - 90)
            x_end_arr = center_x + arrow_len * math.cos(angle_rad_arr)
            y_end_arr = center_y + arrow_len * math.sin(angle_rad_arr)
            canvas.create_line(center_x, center_y, x_end_arr, y_end_arr, fill="red", width=4, arrow=tk.LAST)

            # Legends
            legend_text = (
                f"Departing Azimuth: {alpha1:.2f}° (Blue)\n"
                f"Final Azimuth        : {alpha2:.2f}° (Red)\n"
                f"Total Distance      : {distance_nm:,.2f} NM"
            )
            label = ttk.Label(popup, text=legend_text, font=("Arial", 11, "italic"))
            label.pack(pady=(2, 10))
            popup.transient(parent)
        except Exception as e:
            messagebox.showerror("Graph Error", f"Could not display graph:\n{e}", parent=parent)

    ttk.Button(btn_frame, text="Calculate", command=calculate).grid(row=0, column=0, padx=10)
    ttk.Button(btn_frame, text="Clear", command=clear).grid(row=0, column=1, padx=10)
    ttk.Button(btn_frame, text="Graph", command=show_graph).grid(row=0, column=2, padx=10)

# End of build_gui
