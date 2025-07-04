import tkinter as tk
from tkinter import ttk, messagebox
from rhumb import Rhumb
import math

def build_gui(parent):
    # ========== Title Label ==========
    title = ttk.Label(
        parent,
        text="Enter Coordinates in Degrees and Decimal Minutes",
        font=("Arial", 10, "italic")
    )
    title.pack(pady=(10, 2))
    
    # ========== EXPLANATION ==========
    EXPLANATION_MESSAGE = (
    "\nCalculate the course angle and the distance between the two points on loxodrome (rhumb line)."
    )
    
    # ========== Main Frame ==========
    main_frame = ttk.Frame(parent)
    main_frame.pack(pady=5, padx=5, fill="x", expand=True)

    # ========== Input Section Layout ==========
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

    # ========== Button Frame ==========
    btn_frame = ttk.Frame(parent)
    btn_frame.pack(pady=(5, 0))

    # ========== Result Text ==========
    result_text = tk.Text(
        parent,
        height=12,
        width=66,
        font=("Courier New", 11),
        wrap="word",
        state="normal"  # Start as normal to insert the message
    )
    result_text.insert(tk.END, EXPLANATION_MESSAGE)
    result_text.config(state="disabled")
    result_text.pack(pady=5)

    # ========== Calculation and Utility Functions ==========
    def ddm_to_decimal(degrees, minutes, hemi):
        try:
            degrees = float(degrees)
            minutes = float(minutes)
            value = degrees + minutes / 60
            if hemi in ("S", "W"):
                value = -value
            return value
        except Exception:
            raise ValueError("Degrees and minutes must be numbers.")

    def decimal_to_ddm(val, hemi_pos, hemi_neg, deg_digits=2):
        hemi = hemi_pos if val >= 0 else hemi_neg
        abs_val = abs(val)
        degrees = int(abs_val)
        minutes = (abs_val - degrees) * 60
        deg_str = str(degrees).zfill(deg_digits)
        min_str = f"{minutes:05.2f}"
        return f"{deg_str}° {min_str}' {hemi}"

    def get_decimal(deg_entry, min_entry, hemi_var):
        deg = deg_entry.get()
        mins = min_entry.get()
        hemi = hemi_var.get()
        if not deg or not mins:
            raise ValueError("Both degrees and minutes are required.")
        return ddm_to_decimal(deg, mins, hemi)

    def quadrant_from_azimuth(az):
        if 0 <= az < 90:
            return "NE - 1st Quadrant"
        elif 90 <= az < 180:
            return "SE - 2nd Quadrant"
        elif 180 <= az < 270:
            return "SW - 3rd Quadrant"
        else:
            return "NW - 4th Quadrant"

    latest_results = {"azimuth": None, "distance_nm": None}

    def calculate():
        try:
            # Get start/arrival coordinates in decimal degrees
            lat1 = get_decimal(dep_lat_deg, dep_lat_min, dep_lat_hemi)
            lon1 = get_decimal(dep_lon_deg, dep_lon_min, dep_lon_hemi)
            lat2 = get_decimal(arr_lat_deg, arr_lat_min, arr_lat_hemi)
            lon2 = get_decimal(arr_lon_deg, arr_lon_min, arr_lon_hemi)

            r = Rhumb()
            res = r.Inverse(lat1, lon1, lat2, lon2)
            s12 = res['s12']
            azi12 = res['azi12']

            distance_nm = s12 / 1852.0
            azimuth_rounded = round(azi12, 1)
            distance_rounded = round(distance_nm, 1)
            quadrant = quadrant_from_azimuth(azi12)

            # DDM output
            lat1_ddm = decimal_to_ddm(lat1, "N", "S", deg_digits=2)
            lon1_ddm = decimal_to_ddm(lon1, "E", "W", deg_digits=3)
            lat2_ddm = decimal_to_ddm(lat2, "N", "S", deg_digits=2)
            lon2_ddm = decimal_to_ddm(lon2, "E", "W", deg_digits=3)

            result_str = (
                f"\n"
                f"   --- Heading & Distance Calculator Results ---\n\n"
                f"   Azimuth : {azimuth_rounded:.2f}° ({quadrant})\n"
                f"   Distance: {distance_nm:.2f} NM"                
            )
            result_text.config(state="normal")
            result_text.delete("1.0", tk.END)
            result_text.insert(tk.END, result_str)
            result_text.config(state="disabled")

            # Store for graph
            latest_results["azimuth"] = azi12
            latest_results["distance_nm"] = distance_nm

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
        result_text.config(state="normal")
        result_text.delete("1.0", tk.END)
        result_text.insert(tk.END, EXPLANATION_MESSAGE)   # <-- Insert message again
        result_text.config(state="disabled")
        latest_results["azimuth"] = None
        latest_results["distance_nm"] = None

    def show_graph():
        try:
            azimuth = latest_results["azimuth"]
            distance_nm = latest_results["distance_nm"]
            if azimuth is None:
                messagebox.showinfo("Graph", "Please perform a calculation first.", parent=parent)
                return

            popup = tk.Toplevel(parent)
            popup.title("Heading & Distance Graph")
            try:
                popup.iconbitmap("compass256.ico")
            except Exception:
                pass
            popup.geometry("350x410")
            popup.resizable(False, False)

            canvas = tk.Canvas(popup, width=330, height=330, bg="white")
            canvas.pack(pady=10)

            center_x = 165
            center_y = 165
            radius = 120

            # Draw the earth circle
            canvas.create_oval(center_x - radius, center_y - radius, center_x + radius, center_y + radius, outline="black", width=2)
            font_c = ("Arial", 12, "bold")
            canvas.create_text(center_x, center_y - radius - 15, text="N", font=font_c)
            canvas.create_text(center_x + radius + 15, center_y, text="E", font=font_c)
            canvas.create_text(center_x, center_y + radius + 15, text="S", font=font_c)
            canvas.create_text(center_x - radius - 15, center_y, text="W", font=font_c)
            canvas.create_line(center_x, center_y - radius, center_x, center_y + radius, dash=(2, 2))
            canvas.create_line(center_x - radius, center_y, center_x + radius, center_y, dash=(2, 2))

            # Draw red azimuth arrow
            arrow_len = radius * 0.92
            angle_rad = math.radians(azimuth - 90)
            x_end = center_x + arrow_len * math.cos(angle_rad)
            y_end = center_y + arrow_len * math.sin(angle_rad)
            canvas.create_line(center_x, center_y, x_end, y_end, fill="red", width=4, arrow=tk.LAST)

            # Legend
            legend = (
                f"Azimuth  : {azimuth:.2f}°\n"
                f"Distance: {distance_nm:.2f} NM"
            )
            label = ttk.Label(popup, text=legend, font=("Arial", 11, "italic"))
            label.pack(pady=(2, 10))
            popup.transient(parent)
        except Exception as e:
            messagebox.showerror("Graph Error", f"Could not display graph:\n{e}", parent=parent)

    # ========== Buttons ==========
    ttk.Button(btn_frame, text="Calculate", command=calculate).grid(row=0, column=0, padx=10)
    ttk.Button(btn_frame, text="Clear", command=clear).grid(row=0, column=1, padx=10)
    ttk.Button(btn_frame, text="Graph", command=show_graph).grid(row=0, column=2, padx=10)

# --- End of build_gui ---
