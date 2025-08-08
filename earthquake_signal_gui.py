import dearpygui.dearpygui as dpg
import scipy.io
from scipy import signal as scipy_signal
import numpy as np
import os
from scipy.io import savemat
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from PIL import Image
import math
import matplotlib
import pandas as pd
matplotlib.use('Agg')

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0  # Dünya yarıçapı (km)
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1)*math.cos(phi2)*math.sin(dlambda/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return R * c

def sp_distance_from_ptime_stime(ptime, stime, vp=6.0, vs=3.0):
    sp_time = stime - ptime
    return sp_time * (vp * vs) / (vp - vs)

def calculate_azimuth(lat1, lon1, lat2, lon2):
    import math
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)
    dlon = math.radians(lon2 - lon1)
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    azimuth = math.atan2(x, y)
    azimuth = math.degrees(azimuth)
    azimuth = (azimuth + 360) % 360
    return azimuth

def detect_p_polarity(signal, p_index, window=5):
    import numpy as np
    
    # Extract a larger window around P arrival for detrending
    start_idx = max(0, p_index - window*4)
    end_idx = min(len(signal), p_index + window*4)
    signal_window = signal[start_idx:end_idx]
    
    # Detrend the signal window to remove DC offset and linear trends
    detrended_window = scipy_signal.detrend(signal_window, type='linear')
    
    # Adjust p_index for the windowed signal
    adjusted_p_index = p_index - start_idx
    
    if adjusted_p_index < len(detrended_window) and adjusted_p_index >= window:
        # Calculate pre-event noise baseline from detrended signal
        noise_start = max(0, adjusted_p_index - window*2)  # Further back for noise
        noise_end = max(0, adjusted_p_index - window//2)   # Stop before P arrival
        pre_event_noise = detrended_window[noise_start:noise_end]
        
        if len(pre_event_noise) > 0:
            noise_baseline = np.mean(pre_event_noise)
        else:
            noise_baseline = 0.0
        
        # Get P arrival value from detrended signal
        p_value = detrended_window[adjusted_p_index]
        
        # Calculate difference from pre-event noise baseline
        diff = p_value - noise_baseline
        
        # Simple sign check - no threshold needed after detrend + noise baseline
        if diff > 0:
            return 'up'      # Positive = up motion
        elif diff < 0:
            return 'down'    # Negative = down motion
        else:
            return 'flat'    # Exactly zero (extremely rare)
    else:
        return None

def calculate_signal_azimuth(signals, p_index):
    # North = Channel 1 (0), East = Channel 2 (1)
    # Use the passed signal (can be raw or detrended)
    sig = signals[0] if isinstance(signals, list) else signals
    north_val = sig[p_index, 0]
    east_val = sig[p_index, 1]
    
    # DEBUG: Print values before calculation
    print(f"[DEBUG] calc_sig_azimuth: north_val={north_val:.6f}, east_val={east_val:.6f}")
    
    # Use atan2 for direct 4-quadrant azimuth calculation  
    azimuth_rad = np.arctan2(east_val, north_val)
    azimuth_deg = np.degrees(azimuth_rad)
    if azimuth_deg < 0:
        azimuth_deg += 360
    
    print(f"[DEBUG] calc_sig_azimuth: atan2({east_val:.6f}, {north_val:.6f}) = {azimuth_deg:.2f}°")
    
    return azimuth_deg

def destination_point(lat, lon, distance_km, azimuth_deg):
    import numpy as np
    R = 6371.0  # Dünya yarıçapı (km)
    azimuth_rad = np.radians(azimuth_deg)
    lat1 = np.radians(lat)
    lon1 = np.radians(lon)
    d_div_r = distance_km / R
    lat2 = np.arcsin(np.sin(lat1) * np.cos(d_div_r) +
                     np.cos(lat1) * np.sin(d_div_r) * np.cos(azimuth_rad))
    lon2 = lon1 + np.arctan2(np.sin(azimuth_rad) * np.sin(d_div_r) * np.cos(lat1),
                             np.cos(d_div_r) - np.sin(lat1) * np.sin(lat2))
    return np.degrees(lat2), np.degrees(lon2)

def get_file_status(folder, filename):
    """Get P and S time status for a file and return status indicator"""
    try:
        mat = scipy.io.loadmat(os.path.join(folder, filename), struct_as_record=False, squeeze_me=True)
        
        # Get anEQ structure
        anEQ = None
        if 'EQ' in mat and hasattr(mat['EQ'], 'anEQ'):
            anEQ = mat['EQ'].anEQ
        elif 'anEQ' in mat:
            anEQ = mat['anEQ']
        
        if anEQ is None:
            return "[--]"
        
        # Check P and S times
        ptime_valid = hasattr(anEQ, 'Ptime') and anEQ.Ptime not in [None, -1] and anEQ.Ptime > 0
        stime_valid = hasattr(anEQ, 'Stime') and anEQ.Stime not in [None, -1] and anEQ.Stime > 0
        
        if ptime_valid and stime_valid:
            return "[PS]"
        elif ptime_valid and not stime_valid:
            return "[P-]"
        elif not ptime_valid and stime_valid:
            return "[-S]"
        else:
            return "[--]"
            
    except Exception as e:
        print(f"[DEBUG] Error reading status for {filename}: {e}")
        return "[??]"

def search_files_callback():
    """Search and filter file list based on search term"""
    update_file_list_with_status()

def update_file_list_with_status():
    """Update file list with P/S status indicators"""
    folder = dpg.get_value("Selected Folder")
    files = [f for f in os.listdir(folder) if f.endswith('.mat')]
    
    # Apply search filter if search term exists
    search_term = dpg.get_value("File Search").lower().strip()
    if search_term:
        files = [f for f in files if search_term in f.lower()]
    
    # Remember current selection
    current_selection = dpg.get_value("File List")
    current_filename = None
    if current_selection and " [" in current_selection:
        current_filename = current_selection.split(" [")[0]
    elif current_selection:
        current_filename = current_selection
    
    # Create display list with status indicators
    display_files = []
    selected_index = 0
    for i, filename in enumerate(files):
        status = get_file_status(folder, filename)
        display_name = f"{filename} {status}"
        display_files.append(display_name)
        
        # Find the index of previously selected file
        if current_filename and filename == current_filename:
            selected_index = i
    
    # Update listbox
    dpg.configure_item("File List", items=display_files)
    
    # Update file count with search info
    search_term = dpg.get_value("File Search").lower().strip()
    if search_term:
        dpg.set_value("FileCountText", f"Found {len(files)} files (filtered)")
    else:
        dpg.set_value("FileCountText", f"Files in folder: {len(files)}")
    
    # Restore selection
    if display_files:
        if selected_index < len(display_files):
            dpg.set_value("File List", display_files[selected_index])
        else:
            dpg.set_value("File List", display_files[0])

DEFAULT_FOLDER = "sample_data"
selected_file = [None]
signals = [None]
detrended_signals = [None]  # Store detrended versions for visualization

# State for cursor and picks
cursor_x = [0]
p_pick = [None]
s_pick = [None]

# Helper to update plot titles with file and channel info
def update_plot_titles(filename, sig):
    orientations = ["North", "East", "Up"]
    for i in range(3):
        orientation = orientations[i]
        if sig is not None and i < sig.shape[1]:
            dpg.configure_item(f"PlotWin##{i}", label=f"Channel {i+1} ({orientation}) - {filename}")
        else:
            dpg.configure_item(f"PlotWin##{i}", label=f"Channel {i+1} ({orientation}) - (no data)")

sync_lock = [False]

def sync_cursor_p(sender, app_data):
    tag = dpg.get_item_alias(sender) if isinstance(sender, int) else sender
    if sync_lock[0]:
        return
    sync_lock[0] = True
    x = dpg.get_value(sender)
    p_pick[0] = x
    dpg.set_value("P Pick", f"{x / 100:.2f}")
    
    # Update P cursor Y values for all 3 channels (using detrended signals)
    display_sig = detrended_signals[0] if detrended_signals[0] is not None else signals[0]
    if display_sig is not None:
        try:
            index = int(round(x))
            if 0 <= index < display_sig.shape[0]:
                for ch in range(min(3, display_sig.shape[1])):
                    y_value = display_sig[index, ch]
                    dpg.set_value(f"P Y{ch}", f"{y_value:.6f}")

                
                # Fill remaining channels if signal has fewer than 3 channels
                for ch in range(display_sig.shape[1], 3):
                    dpg.set_value(f"P Y{ch}", "N/A")
            else:
                for ch in range(3):
                    dpg.set_value(f"P Y{ch}", "N/A")
        except:
            for ch in range(3):
                dpg.set_value(f"P Y{ch}", "N/A")
    else:
        for ch in range(3):
            dpg.set_value(f"P Y{ch}", "N/A")
    
    for i in range(3):
        tag2 = f"PPickLine##{i}"
        if tag2 != tag:
            dpg.set_value(tag2, x)
    
    # Update polarity and azimuth when P cursor moves
    update_azimuth_calc_boxes()
    
    sync_lock[0] = False

def sync_cursor_s(sender, app_data):
    tag = dpg.get_item_alias(sender) if isinstance(sender, int) else sender
    if sync_lock[0]:
        return
    sync_lock[0] = True
    x = dpg.get_value(sender)
    s_pick[0] = x
    dpg.set_value("S Pick", f"{x / 100:.2f}")
    
    # Update S cursor Y values for all 3 channels (using detrended signals)
    display_sig = detrended_signals[0] if detrended_signals[0] is not None else signals[0]
    if display_sig is not None:
        try:
            index = int(round(x))
            if 0 <= index < display_sig.shape[0]:
                for ch in range(min(3, display_sig.shape[1])):
                    y_value = display_sig[index, ch]
                    dpg.set_value(f"S Y{ch}", f"{y_value:.6f}")
                # Fill remaining channels if signal has fewer than 3 channels
                for ch in range(display_sig.shape[1], 3):
                    dpg.set_value(f"S Y{ch}", "N/A")
            else:
                for ch in range(3):
                    dpg.set_value(f"S Y{ch}", "N/A")
        except:
            for ch in range(3):
                dpg.set_value(f"S Y{ch}", "N/A")
    else:
        for ch in range(3):
            dpg.set_value(f"S Y{ch}", "N/A")
    
    for i in range(3):
        tag2 = f"SPickLine##{i}"
        if tag2 != tag:
            dpg.set_value(tag2, x)
    
    # Update azimuth when S cursor moves (S-P distance changes)
    update_azimuth_calc_boxes()
    
    sync_lock[0] = False

# --- Cartopy map drawing ---
def draw_turkey_map(epicenter=None, statco=None, out_file="turkey_map.png", azimuth_calc=None, sp_distance=None, event_azimuth_point=None):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import os

    print(f"[DEBUG] draw_turkey_map called with epicenter={epicenter}, statco={statco}, azimuth_calc={azimuth_calc}, sp_distance={sp_distance}, event_azimuth_point={event_azimuth_point}, out_file={os.path.abspath(out_file)}")
    fig = plt.figure(figsize=(5.5, 3.5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([25, 46, 35, 43])
    ax.add_feature(cfeature.BORDERS, linewidth=1)
    ax.add_feature(cfeature.COASTLINE, linewidth=1)
    ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.3, linestyle='--')

    # Epicenter: küçük mavi dolu daire
    if epicenter is not None and len(epicenter) == 2:
        print(f"[DEBUG] Plotting epicenter at lon={epicenter[1]}, lat={epicenter[0]}")
        ax.scatter([epicenter[1]], [epicenter[0]], s=50, color='blue', marker='o', zorder=30, transform=ccrs.PlateCarree())
    # Statco: küçük kırmızı üçgen
    if statco is not None and len(statco) == 2:
        print(f"[DEBUG] Plotting statco at lon={statco[1]}, lat={statco[0]}")
        ax.scatter([statco[1]], [statco[0]], s=50, color='red', marker='v', zorder=31, transform=ccrs.PlateCarree())
    # Predicted epicenter: yeşil daire
    import math
    print(f"[DEBUG] statco={statco} ({type(statco)}), azimuth_calc={azimuth_calc} ({type(azimuth_calc)}), sp_distance={sp_distance} ({type(sp_distance)})")
    valid = (
        statco is not None and azimuth_calc is not None and sp_distance is not None and
        isinstance(statco, (list, tuple, np.ndarray)) and len(statco) == 2 and
        all([isinstance(x, (int, float, np.floating, np.integer)) and not math.isnan(x) and not math.isinf(x) for x in statco]) and
        isinstance(azimuth_calc, (int, float, np.floating, np.integer)) and not math.isnan(azimuth_calc) and not math.isinf(azimuth_calc) and
        isinstance(sp_distance, (int, float, np.floating, np.integer)) and not math.isnan(sp_distance) and not math.isinf(sp_distance)
    )
    if valid:
        pred_lat, pred_lon = destination_point(statco[0], statco[1], sp_distance, azimuth_calc)
        print(f"[DEBUG] Plotting predicted epicenter at lon={pred_lon}, lat={pred_lat}")
        ax.scatter([pred_lon], [pred_lat], s=50, color='green', marker='o', zorder=100, transform=ccrs.PlateCarree())
    else:
        print("[DEBUG] Predicted epicenter çizilmiyor: Değerlerden biri geçersiz veya yanlış tipte.")
    # --- Legend: sağ alt köşe ---
    legend_x = 37.8  # longitude (biraz daha sola kaydırıldı)
    legend_y = 36.3  # latitude (biraz yukarı taşıdık)
    ax.scatter([legend_x], [legend_y], s=40, color='blue', marker='o', zorder=40, transform=ccrs.PlateCarree())
    ax.text(legend_x + 0.3, legend_y, "Epicenter", color='black', fontsize=10, va='center', zorder=41, transform=ccrs.PlateCarree())
    ax.scatter([legend_x], [legend_y-0.5], s=40, color='red', marker='v', zorder=40, transform=ccrs.PlateCarree())
    ax.text(legend_x + 0.3, legend_y-0.5, "Station", color='black', fontsize=10, va='center', zorder=41, transform=ccrs.PlateCarree())
    ax.scatter([legend_x], [legend_y-1.0], s=40, color='green', marker='o', zorder=40, transform=ccrs.PlateCarree())
    ax.text(legend_x + 0.3, legend_y-1.0, "Predicted Epicenter", color='black', fontsize=10, va='center', zorder=41, transform=ccrs.PlateCarree())

    plt.tight_layout()
    plt.savefig(out_file, dpi=100)
    plt.close(fig)
    import time
    time.sleep(0.2)  # PNG'nin tam yazılması için bekleme süresi artırıldı
    print(f"[DEBUG] draw_turkey_map: PNG saved at {out_file}")
    print(f"[DEBUG] Map PNG updated: {out_file}")

# initialize_file_list fonksiyonunu eski haline döndür:
def initialize_file_list():
    folder = DEFAULT_FOLDER
    dpg.set_value("Selected Folder", folder)
    dpg.set_value("File Search", "")  # Clear search box
    update_file_list_with_status()
    
    # Listbox'taki gerçek item sayısını kontrol et
    listbox_items = dpg.get_item_configuration("File List")["items"]
    print(f"[DEBUG] Listbox contains {len(listbox_items)} items")
    
    if listbox_items:
        dpg.set_value("File List", listbox_items[0])

# Listbox render'ında seçili dosyanın rengini güncelle
# (DearPyGui Listbox doğrudan renkli item desteklemez, seçili dosya adını üstte renkle gösterebiliriz)
# Dosya okuma ve veri erişiminde hem mat['EQ'].anEQ.Accel hem de mat['anEQ'].Accel desteklensin:
def file_select_callback(sender, app_data):
    import numpy as np
    folder = dpg.get_value("Selected Folder")
    display_filename = dpg.get_value("File List")
    
    # Extract actual filename from display (remove status indicator)
    if " [" in display_filename:
        filename = display_filename.split(" [")[0]
    else:
        filename = display_filename
        
    print(f"[DEBUG] file_select_callback: Dosya seçildi: {filename}")
    
    # Seçili dosyanın listede kaçıncı sırada olduğunu bul
    all_display_files = dpg.get_item_configuration("File List")["items"]
    if display_filename in all_display_files:
        file_index = all_display_files.index(display_filename) + 1
        dpg.set_value("FileCountText", f"File {file_index} of {len(all_display_files)}: {filename}")
    
    try:
        mat = scipy.io.loadmat(os.path.join(folder, filename), struct_as_record=False, squeeze_me=True)
        # Accel erişimi için iki farklı yol dene
        sig = None
        anEQ = None
        if 'EQ' in mat and hasattr(mat['EQ'], 'anEQ') and hasattr(mat['EQ'].anEQ, 'Accel'):
            sig = mat['EQ'].anEQ.Accel
            anEQ = mat['EQ'].anEQ
            print('[DEBUG] Veri EQ.anEQ.Accel üzerinden okundu')
        elif 'anEQ' in mat and hasattr(mat['anEQ'], 'Accel'):
            sig = mat['anEQ'].Accel
            anEQ = mat['anEQ']
            print('[DEBUG] Veri anEQ.Accel üzerinden okundu')
        else:
            print('[DEBUG] Uygun Accel alanı bulunamadı!')
            sig = None
            anEQ = None
        signals[0] = sig
        
        # Calculate detrended version for visualization (doesn't modify original data)
        if sig is not None and isinstance(sig, np.ndarray) and sig.ndim >= 2:
            detrended_sig = np.zeros_like(sig)
            for ch in range(sig.shape[1]):
                detrended_sig[:, ch] = scipy_signal.detrend(sig[:, ch], type='linear')
            detrended_signals[0] = detrended_sig
        else:
            detrended_signals[0] = None
        
        selected_file[0] = os.path.join(folder, filename)
        # Reset picks and drag lines
        # Ptime ve Stime varsa, cursorları oradan başlat
        p_val = float(anEQ.Ptime) * 100 if hasattr(anEQ, 'Ptime') and anEQ.Ptime not in [None, -1] else 0.0
        s_val = float(anEQ.Stime) * 100 if hasattr(anEQ, 'Stime') and anEQ.Stime not in [None, -1] else 0.0
        p_pick[0] = p_val
        s_pick[0] = s_val
        dpg.set_value("P Pick", f"{p_val / 100:.2f}")
        dpg.set_value("S Pick", f"{s_val / 100:.2f}")
        
        # Initialize Y values for P and S cursors using detrended signals
        display_sig = detrended_signals[0] if detrended_signals[0] is not None else sig
        if isinstance(display_sig, np.ndarray) and display_sig.ndim >= 2:
            try:
                p_index = int(round(p_val))
                s_index = int(round(s_val))
                
                # Initialize P Y values for all channels (detrended)
                if 0 <= p_index < display_sig.shape[0]:
                    for ch in range(min(3, display_sig.shape[1])):
                        p_y_value = display_sig[p_index, ch]
                        dpg.set_value(f"P Y{ch}", f"{p_y_value:.6f}")
                    # Fill remaining channels if signal has fewer than 3 channels
                    for ch in range(display_sig.shape[1], 3):
                        dpg.set_value(f"P Y{ch}", "N/A")
                else:
                    for ch in range(3):
                        dpg.set_value(f"P Y{ch}", "N/A")
                        
                # Initialize S Y values for all channels (detrended)
                if 0 <= s_index < display_sig.shape[0]:
                    for ch in range(min(3, display_sig.shape[1])):
                        s_y_value = display_sig[s_index, ch]
                        dpg.set_value(f"S Y{ch}", f"{s_y_value:.6f}")
                    # Fill remaining channels if signal has fewer than 3 channels
                    for ch in range(display_sig.shape[1], 3):
                        dpg.set_value(f"S Y{ch}", "N/A")
                else:
                    for ch in range(3):
                        dpg.set_value(f"S Y{ch}", "N/A")
            except:
                for ch in range(3):
                    dpg.set_value(f"P Y{ch}", "N/A")
                    dpg.set_value(f"S Y{ch}", "N/A")
        else:
            for ch in range(3):
                dpg.set_value(f"P Y{ch}", "N/A")
                dpg.set_value(f"S Y{ch}", "N/A")
        
        # Use same display signal for both plots and Y boxes consistency
        display_sig_for_plot = detrended_signals[0] if detrended_signals[0] is not None else sig
        
        for i in range(3):
            # Use the same signal as Y boxes for consistency
            display_sig = display_sig_for_plot
            if isinstance(display_sig, np.ndarray) and display_sig.ndim >= 2 and i < display_sig.shape[1]:
                x = np.arange(display_sig.shape[0])
                y = display_sig[:, i]
                dpg.set_value(f"Plot##{i}", (x.tolist(), y.tolist()))
                dpg.set_value(f"PPickLine##{i}", p_val)
                dpg.set_value(f"SPickLine##{i}", s_val)
            else:
                dpg.set_value(f"Plot##{i}", ([], []))
                dpg.set_value(f"PPickLine##{i}", p_val)
                dpg.set_value(f"SPickLine##{i}", s_val)
        update_plot_titles(filename, sig)
        # Only show selected fields in anEQ, no type/shape
        try:
            # anEQ yukarıda seçildiği gibi kullanılacak
            show_fields = [
                'place', 'date', 'epicenter', 'depth', 'magnitudeval', 'statID', 'magnitudetype',
                'statco', 'numofData', 'pga', 'Vs30', 'Ptime', 'Stime'
            ]
            attrs = []
            if anEQ is not None:
                for attr in show_fields:
                    if hasattr(anEQ, attr):
                        val = getattr(anEQ, attr)
                        preview = ""
                        if isinstance(val, np.ndarray):
                            if val.size == 0:
                                preview = "[]"
                            elif val.ndim == 0:
                                preview = f"{val.item()}"
                            elif val.ndim == 1:
                                preview = np.array2string(val[:5], separator=', ', precision=3, max_line_width=60)
                                if val.size > 5:
                                    preview = preview.rstrip(']') + ', ...]'
                            elif val.ndim == 2:
                                rows = min(5, val.shape[0])
                                preview = np.array2string(val[:rows], separator=', ', precision=3, max_line_width=60)
                                if val.shape[0] > 5:
                                    preview = preview.rstrip(']') + ', ...]'
                            else:
                                preview = f"array shape {val.shape}"
                        else:
                            try:
                                # Special formatting for Ptime and Stime - show only 2 decimal places
                                if attr in ['Ptime', 'Stime'] and isinstance(val, (int, float, np.floating, np.integer)):
                                    preview = f"{float(val):.2f}"
                                else:
                                    preview = str(val)
                            except Exception:
                                preview = "?"
                        attrs.append(f"{attr}: {preview}")
                dpg.set_value("anEQ Attributes List", "\n".join(attrs) if attrs else "No selected attributes found.")
            else:
                dpg.set_value("anEQ Attributes List", "No anEQ structure found in MAT file.")
        except Exception as e:
            dpg.set_value("anEQ Attributes List", f"Error reading attributes: {e}")

        # --- Extract epicenter/statco and update distance box ---
        epicenter = None
        statco = None
        ptime = None
        stime = None
        azimuth = None
        back_azimuth = None
        p_polarity = None
        event_azimuth = None
        if hasattr(anEQ, 'epicenter') and isinstance(anEQ.epicenter, (np.ndarray, list, tuple)) and len(anEQ.epicenter) == 2:
            epicenter = [float(anEQ.epicenter[0]), float(anEQ.epicenter[1])]
        if hasattr(anEQ, 'statco') and isinstance(anEQ.statco, (np.ndarray, list, tuple)) and len(anEQ.statco) == 2:
            statco = [float(anEQ.statco[0]), float(anEQ.statco[1])]
        if hasattr(anEQ, 'Ptime'):
            try:
                ptime = float(anEQ.Ptime)
            except Exception:
                ptime = None
        if hasattr(anEQ, 'Stime'):
            try:
                stime = float(anEQ.Stime)
            except Exception:
                stime = None
        print(f"[DEBUG] epicenter: {epicenter}, statco: {statco}, len(epicenter): {len(epicenter) if epicenter is not None else 'None'}, len(statco): {len(statco) if statco is not None else 'None'}")
        # Orijinal azimuth sadece coğrafi olarak güncellensin:
        if epicenter is not None and statco is not None and len(epicenter) == 2 and len(statco) == 2:
            try:
                lat1, lon1 = float(epicenter[0]), float(epicenter[1])
                lat2, lon2 = float(statco[0]), float(statco[1])
                print(f"[DEBUG] Calculating distance between (lat1={lat1}, lon1={lon1}) and (lat2={lat2}, lon2={lon2})")
                distance_km = haversine(lat1, lon1, lat2, lon2)
                print(f"[DEBUG] Calculated distance: {distance_km} km")
                dpg.set_value("Epicenter-Statco Distance", f"{distance_km:.2f} km")
                # Azimuth and back azimuth
                azimuth = calculate_azimuth(lat2, lon2, lat1, lon1)  # statco -> epicenter
                back_azimuth = calculate_azimuth(lat1, lon1, lat2, lon2)  # epicenter -> statco
                print(f"[DEBUG] Calculated azimuth: {azimuth}")
                dpg.set_value("Azimuth", f"{azimuth:.2f}°")  # SADECE burada güncelleniyor
                dpg.set_value("Back Azimuth", f"{back_azimuth:.2f}°")
            except Exception as e:
                print(f"[DEBUG] Distance calculation error: {e}")
                dpg.set_value("Epicenter-Statco Distance", "-")
                dpg.set_value("Azimuth", "-")
                dpg.set_value("Back Azimuth", "-")
        else:
            dpg.set_value("Epicenter-Statco Distance", "-")
            dpg.set_value("Azimuth", "-")
            dpg.set_value("Back Azimuth", "-")

        # S-P distance calculation (Vp/Vs GUI'den alınır)
        if ptime is not None and stime is not None and ptime > 0 and stime > 0 and ptime != -1 and stime != -1:
            try:
                vp = float(dpg.get_value("Vp Value"))
                vs = float(dpg.get_value("Vs Value"))
                sp_dist = sp_distance_from_ptime_stime(ptime, stime, vp=vp, vs=vs)
                dpg.set_value("SP Distance", f"{sp_dist:.2f} km")
            except Exception as e:
                print(f"[DEBUG] S-P distance calculation error: {e}")
                dpg.set_value("SP Distance", "-")
        else:
            dpg.set_value("SP Distance", "-")

        # P polarity and event azimuth
        if ptime is not None and ptime > 0 and ptime != -1 and signals[0] is not None:
            try:
                # Use channel 2 (Up) for polarity - dikey hareket daha net
                signal = signals[0][:, 2]
                sample_rate = 100  # Varsayım, gerekirse değiştir
                p_index = int(ptime * sample_rate)
                p_polarity = detect_p_polarity(signal, p_index, window=5)
                dpg.set_value("P Polarity", p_polarity if p_polarity else "-")
                # Event azimuth: sadece coğrafi azimuth (statco->epicenter) yazılsın, polarity'ye göre değişmesin
                if azimuth is not None:
                    event_azimuth = azimuth
                else:
                    event_azimuth = None
                dpg.set_value("Event Azimuth", f"{event_azimuth:.2f}°" if event_azimuth is not None else "-")
            except Exception as e:
                print(f"[DEBUG] P polarity/azimuth error: {e}")
                dpg.set_value("P Polarity", "-")
                dpg.set_value("Event Azimuth", "-")
        else:
            dpg.set_value("P Polarity", "-")
            dpg.set_value("Event Azimuth", "-")

        # Azimuth kutusunu Event Azimuth (Calc.) ile güncelle
        event_az_calc_str = dpg.get_value("Event Azimuth (Calc.)")
        print(f"[DEBUG] Event Azimuth (Calc.) kutusundan okunan: '{event_az_calc_str}'")
        if event_az_calc_str and event_az_calc_str != "-":
            azimuth_calc = float(event_az_calc_str.replace("°", "").strip())
        else:
            azimuth_calc = None

        # Azimuth (Calc.) ve SP Distance değerlerini oku
        try:
            sp_distance = None
            event_az_calc_str = dpg.get_value("Event Azimuth (Calc.)")
            sp_dist_str = dpg.get_value("SP Distance")
            if event_az_calc_str and event_az_calc_str != "-":
                azimuth_calc = float(event_az_calc_str.replace("°", "").strip())
            if sp_dist_str and sp_dist_str != "-":
                sp_distance = float(sp_dist_str.replace("km", "").strip())
        except Exception:
            azimuth_calc = None
            sp_distance = None

        # --- Kutuları güncelle ---
        update_azimuth_calc_boxes()
        # --- Haritayı çiz ---
        try:
            print(f"[DEBUG] Drawing map with epicenter: {epicenter}, statco: {statco}, azimuth_calc: {azimuth_calc}, sp_distance: {sp_distance}")
            # Event azimuth referans noktası için hesaplama
            event_azimuth_point = None
            if statco is not None and azimuth_calc is not None and sp_distance is not None:
                event_azimuth_point = destination_point(statco[0], statco[1], sp_distance, azimuth_calc)
            draw_turkey_map(epicenter, statco, out_file=MAP_PNG_PATH, azimuth_calc=azimuth_calc, sp_distance=sp_distance, event_azimuth_point=event_azimuth_point)
            print(f"[DEBUG] After draw_turkey_map, PNG path: {MAP_PNG_PATH}")
            import time
            time.sleep(0.1)  # PNG'nin tam yazılması için kısa bekleme
            print(f"[DEBUG] Slept 0.1s before loading PNG for texture.")
            from PIL import Image
            print(f"[DEBUG] Loading PNG for GUI: {MAP_PNG_PATH}")
            image = Image.open(MAP_PNG_PATH).convert("RGBA")
            map_data = np.array(image)
            print(f"[DEBUG] PNG shape: {map_data.shape}, dtype: {map_data.dtype}, min: {map_data.min()}, max: {map_data.max()}")
            if map_data.shape[2] != 4:
                raise ValueError("PNG is not RGBA! Shape: {}".format(map_data.shape))
            map_data = map_data.astype(np.float32) / 255.0
            map_data = map_data.copy()  # contiguous
            width, height = image.size  # width, height
            if dpg.does_item_exist("TurkeyMapTex"):
                pass
            print("[DEBUG] add_static_texture called (set_value mode)")
            dpg.set_value("TurkeyMapTex", map_data.flatten().tolist())
            print("[DEBUG] configure_item called for TurkeyMapImage")
            dpg.configure_item("TurkeyMapImage", texture_tag="TurkeyMapTex", width=width, height=height)
        except Exception as e:
            print(f"[DEBUG] Map PNG update error: {e}")

    except Exception as e:
        print(f"[DEBUG] File could not be read or is not in expected format: {e}")
        for i in range(3):
            dpg.set_value(f"Plot##{i}", ([], []))
            dpg.configure_item(f"PlotWin##{i}", label=f"Channel {i+1} - (no data)")
            dpg.set_value(f"PPickLine##{i}", 0.0)
            dpg.set_value(f"SPickLine##{i}", 0.0)
        dpg.set_value("anEQ Attributes List", "-")
        try:
            dpg.delete_item("TurkeyMapTex")
        except: pass

    # P pick ve polarite kutuları güncellendikten sonra:
    update_azimuth_calc_boxes()

def save_excel_callback():
    """Process all files and save data to Excel"""
    try:
        folder = dpg.get_value("Selected Folder")
        files = [f for f in os.listdir(folder) if f.endswith('.mat')]
        
        if not files:
            print("[DEBUG] No .mat files found for Excel export")
            dpg.set_value("Progress Text", "No .mat files found")
            return
            
        print(f"[DEBUG] Processing {len(files)} files for Excel export...")
        
        # Reset progress bar
        dpg.set_value("Excel Progress", 0.0)
        dpg.set_value("Progress Text", f"Starting... (0/{len(files)})")
        dpg.split_frame()  # Update GUI
        
        # Data lists for DataFrame
        data = {
            'Filename': [],
            'P_time': [],
            'S_time': [],
            'Original_Distance_km': [],
            'Event_Azimuth_deg': [],
            'Predicted_Epicenter_Lat': [],
            'Predicted_Epicenter_Lon': [],
            'Station_Lat': [],
            'Station_Lon': [],
            'Original_Epicenter_Lat': [],
            'Original_Epicenter_Lon': []
        }
        
        # Get Vp/Vs values from GUI
        vp = float(dpg.get_value("Vp Value"))
        vs = float(dpg.get_value("Vs Value"))
        
        for i, filename in enumerate(files):
            # Update progress
            progress = (i) / len(files)
            dpg.set_value("Excel Progress", progress)
            dpg.set_value("Progress Text", f"Processing: {filename} ({i}/{len(files)})")
            dpg.split_frame()  # Update GUI immediately
            
            print(f"[DEBUG] Processing file: {filename}")
            
            try:
                # Load mat file
                mat = scipy.io.loadmat(os.path.join(folder, filename), struct_as_record=False, squeeze_me=True)
                
                # Get anEQ structure
                anEQ = None
                sig = None
                if 'EQ' in mat and hasattr(mat['EQ'], 'anEQ') and hasattr(mat['EQ'].anEQ, 'Accel'):
                    anEQ = mat['EQ'].anEQ
                    sig = mat['EQ'].anEQ.Accel
                elif 'anEQ' in mat and hasattr(mat['anEQ'], 'Accel'):
                    anEQ = mat['anEQ']
                    sig = mat['anEQ'].Accel
                
                # Initialize default values
                ptime = None
                stime = None
                distance_km = None
                event_azimuth = None
                pred_lat = None
                pred_lon = None
                station_lat = None
                station_lon = None
                orig_epicenter_lat = None
                orig_epicenter_lon = None
                
                if anEQ is not None:
                    # Get P and S times
                    if hasattr(anEQ, 'Ptime') and anEQ.Ptime not in [None, -1]:
                        ptime = float(anEQ.Ptime)
                    if hasattr(anEQ, 'Stime') and anEQ.Stime not in [None, -1]:
                        stime = float(anEQ.Stime)
                    
                    # Get epicenter and station coordinates
                    epicenter = None
                    statco = None
                    if hasattr(anEQ, 'epicenter') and isinstance(anEQ.epicenter, (np.ndarray, list, tuple)) and len(anEQ.epicenter) == 2:
                        epicenter = [float(anEQ.epicenter[0]), float(anEQ.epicenter[1])]
                        orig_epicenter_lat, orig_epicenter_lon = epicenter
                    if hasattr(anEQ, 'statco') and isinstance(anEQ.statco, (np.ndarray, list, tuple)) and len(anEQ.statco) == 2:
                        statco = [float(anEQ.statco[0]), float(anEQ.statco[1])]
                        station_lat, station_lon = statco
                    
                    # Calculate distance
                    if epicenter is not None and statco is not None:
                        distance_km = haversine(epicenter[0], epicenter[1], statco[0], statco[1])
                    
                    # Calculate event azimuth and predicted epicenter
                    if ptime is not None and stime is not None and ptime > 0 and stime > 0 and sig is not None and statco is not None:
                        try:
                            # Calculate S-P distance
                            sp_distance = sp_distance_from_ptime_stime(ptime, stime, vp=vp, vs=vs)
                            
                            # Calculate signal azimuth at P pick
                            p_index = int(ptime * 100)  # Assuming 100 Hz sampling rate
                            if p_index < sig.shape[0]:
                                # Calculate azimuth from signal using atan2
                                north_val = sig[p_index, 0] if sig.shape[1] > 0 else 0
                                east_val = sig[p_index, 1] if sig.shape[1] > 1 else 0
                                
                                # Direct atan2 calculation
                                azimuth_rad = np.arctan2(east_val, north_val)
                                azimuth_deg = np.degrees(azimuth_rad)
                                if azimuth_deg < 0:
                                    azimuth_deg += 360
                                
                                # Apply polarity correction
                                if sig.shape[1] > 2:  # Up channel available
                                    signal_up = sig[:, 2]
                                    p_polarity = detect_p_polarity(signal_up, p_index, window=5)
                                    
                                    # atan2 gives direct azimuth, then apply polarity correction:
                                    # Down motion (compression) = direct azimuth toward epicenter
                                    # Up motion (dilation) = azimuth + 180° (opposite direction)
                                    if p_polarity == "down":
                                        event_azimuth = azimuth_deg  # Compression: direct azimuth
                                    elif p_polarity == "up":
                                        event_azimuth = (azimuth_deg + 180) % 360  # Dilation: opposite direction
                                    else:
                                        event_azimuth = azimuth_deg  # Default if polarity unclear
                                else:
                                    event_azimuth = azimuth_deg
                                
                                # Calculate predicted epicenter
                                pred_lat, pred_lon = destination_point(statco[0], statco[1], sp_distance, event_azimuth)
                                
                        except Exception as e:
                            print(f"[DEBUG] Error calculating azimuth/predicted epicenter for {filename}: {e}")
                
                # Add data to lists
                data['Filename'].append(filename)
                data['P_time'].append(ptime)
                data['S_time'].append(stime)
                data['Original_Distance_km'].append(distance_km)
                data['Event_Azimuth_deg'].append(event_azimuth)
                data['Predicted_Epicenter_Lat'].append(pred_lat)
                data['Predicted_Epicenter_Lon'].append(pred_lon)
                data['Station_Lat'].append(station_lat)
                data['Station_Lon'].append(station_lon)
                data['Original_Epicenter_Lat'].append(orig_epicenter_lat)
                data['Original_Epicenter_Lon'].append(orig_epicenter_lon)
                
            except Exception as e:
                print(f"[DEBUG] Error processing file {filename}: {e}")
                # Add empty row for failed files
                data['Filename'].append(filename)
                for key in data.keys():
                    if key != 'Filename':
                        data[key].append(None)
        
        # Update progress to show saving
        dpg.set_value("Excel Progress", 1.0)
        dpg.set_value("Progress Text", f"Saving Excel file... ({len(files)}/{len(files)})")
        dpg.split_frame()
        
        # Create DataFrame and save to Excel
        df = pd.DataFrame(data)
        
        # Save to Excel file in the current code directory
        excel_path = os.path.join(os.getcwd(), "earthquake_analysis_results.xlsx")
        df.to_excel(excel_path, index=False, engine='openpyxl')
        
        print(f"[DEBUG] Excel file saved successfully: {excel_path}")
        print(f"[DEBUG] Processed {len(files)} files")
        
        # Final progress update
        dpg.set_value("Progress Text", f"Complete! {len(files)} files processed")
        
    except Exception as e:
        print(f"[DEBUG] Error in save_excel_callback: {e}")
        dpg.set_value("Progress Text", f"Excel export failed: {e}")
        dpg.set_value("Excel Progress", 0.0)

def mark_p_pick_callback():
    x = dpg.get_value("PPickLine##0")
    p_pick[0] = x
    for i in range(3):
        dpg.set_value(f"PPickLine##{i}", x)
    dpg.set_value("P Pick", f"{x / 100:.2f}")
    
    # Update P cursor Y values for all 3 channels
    if signals[0] is not None:
        try:
            index = int(x)
            if 0 <= index < signals[0].shape[0]:
                for ch in range(min(3, signals[0].shape[1])):
                    y_value = signals[0][index, ch]
                    dpg.set_value(f"P Y{ch}", f"{y_value:.6f}")
                # Fill remaining channels if signal has fewer than 3 channels
                for ch in range(signals[0].shape[1], 3):
                    dpg.set_value(f"P Y{ch}", "N/A")
            else:
                for ch in range(3):
                    dpg.set_value(f"P Y{ch}", "N/A")
        except:
            for ch in range(3):
                dpg.set_value(f"P Y{ch}", "N/A")
    else:
        for ch in range(3):
            dpg.set_value(f"P Y{ch}", "N/A")
    
    update_azimuth_calc_boxes()
    file_select_callback(None, None)

def mark_s_pick_callback():
    x = dpg.get_value("SPickLine##0")
    s_pick[0] = x
    for i in range(3):
        dpg.set_value(f"SPickLine##{i}", x)
    dpg.set_value("S Pick", f"{x / 100:.2f}")
    
    # Update S cursor Y values for all 3 channels
    if signals[0] is not None:
        try:
            index = int(x)
            if 0 <= index < signals[0].shape[0]:
                for ch in range(min(3, signals[0].shape[1])):
                    y_value = signals[0][index, ch]
                    dpg.set_value(f"S Y{ch}", f"{y_value:.6f}")
                # Fill remaining channels if signal has fewer than 3 channels
                for ch in range(signals[0].shape[1], 3):
                    dpg.set_value(f"S Y{ch}", "N/A")
            else:
                for ch in range(3):
                    dpg.set_value(f"S Y{ch}", "N/A")
        except:
            for ch in range(3):
                dpg.set_value(f"S Y{ch}", "N/A")
    else:
        for ch in range(3):
            dpg.set_value(f"S Y{ch}", "N/A")
    
    update_azimuth_calc_boxes()
    file_select_callback(None, None)

def save_ptime_callback():
    path = selected_file[0]
    if not path:
        print("[DEBUG] No file selected for saving Ptime.")
        return
    try:
        dpg.split_frame()  # GUI eventlerini işle, en güncel P Pick değerini al
        p_sec = float(dpg.get_value("P Pick"))
        mat = scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)
        yazildi = False
        if 'EQ' in mat and hasattr(mat['EQ'], 'anEQ') and hasattr(mat['EQ'].anEQ, '__dict__'):
            mat['EQ'].anEQ.Ptime = p_sec
            yazildi = True
        if 'anEQ' in mat and hasattr(mat['anEQ'], '__dict__'):
            mat['anEQ'].Ptime = p_sec
            yazildi = True
        if not yazildi:
            print("[DEBUG] Could not find anEQ in mat file.")
            return
        savemat(path, mat, do_compression=True)
        print(f"[DEBUG] Saved Ptime={p_sec} to {path}")
        
        # Update file list with new status
        update_file_list_with_status()
        file_select_callback(None, None)
    except Exception as e:
        print(f"[DEBUG] Error saving Ptime: {e}")

def save_stime_callback():
    path = selected_file[0]
    if not path:
        print("[DEBUG] No file selected for saving Stime.")
        return
    try:
        dpg.split_frame()  # GUI eventlerini işle, en güncel S Pick değerini al
        s_sec = float(dpg.get_value("S Pick"))
        mat = scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)
        yazildi = False
        if 'EQ' in mat and hasattr(mat['EQ'], 'anEQ') and hasattr(mat['EQ'].anEQ, '__dict__'):
            mat['EQ'].anEQ.Stime = s_sec
            yazildi = True
        if 'anEQ' in mat and hasattr(mat['anEQ'], '__dict__'):
            mat['anEQ'].Stime = s_sec
            yazildi = True
        if not yazildi:
            print("[DEBUG] Could not find anEQ in mat file.")
            return
        savemat(path, mat, do_compression=True)
        print(f"[DEBUG] Saved Stime={s_sec} to {path}")
        
        # Update file list with new status
        update_file_list_with_status()
        file_select_callback(None, None)
    except Exception as e:
        print(f"[DEBUG] Error saving Stime: {e}")

def sync_all_lines():
    # P lines
    p_val = dpg.get_value("PPickLine##0")
    for i in range(1, 3):
        if dpg.get_value(f"PPickLine##{i}") != p_val:
            dpg.set_value(f"PPickLine##{i}", p_val)
    dpg.set_value("P Pick", f"{p_val / 100:.2f}")
    # S lines
    s_val = dpg.get_value("SPickLine##0")
    for i in range(1, 3):
        if dpg.get_value(f"SPickLine##{i}") != s_val:
            dpg.set_value(f"SPickLine##{i}", s_val)
    dpg.set_value("S Pick", f"{s_val / 100:.2f}")

def update_azimuth_calc_boxes():
    try:
        p_sec = float(dpg.get_value("P Pick"))
        p_index = int(p_sec * 100)
        # Use detrended signals for azimuth calculation consistency
        calc_sig = detrended_signals[0] if detrended_signals[0] is not None else signals[0]
        
        # DEBUG: Check which signal is being used
        if detrended_signals[0] is not None:
            print(f"[DEBUG] Using DETRENDED signal for azimuth calc")
        else:
            print(f"[DEBUG] Using RAW signal for azimuth calc (detrended_signals[0] is None)")
        
        if calc_sig is not None and p_index < calc_sig.shape[0]:
            # DEBUG: Check values in calc_sig at p_index
            debug_north = calc_sig[p_index, 0] if calc_sig.shape[1] > 0 else "N/A"
            debug_east = calc_sig[p_index, 1] if calc_sig.shape[1] > 1 else "N/A"
            print(f"[DEBUG] calc_sig values at p_index {p_index}: north={debug_north:.6f}, east={debug_east:.6f}")
            
            # DEBUG: Compare with GUI box values  
            gui_north = dpg.get_value("P Y0")
            gui_east = dpg.get_value("P Y1")
            print(f"[DEBUG] GUI box values: north={gui_north}, east={gui_east}")
            
            # Use GUI box values directly for consistency (these are guaranteed to be detrended)
            gui_north_str = dpg.get_value("P Y0")
            gui_east_str = dpg.get_value("P Y1")
            
            try:
                gui_north = float(gui_north_str) if gui_north_str != "N/A" else 0.0
                gui_east = float(gui_east_str) if gui_east_str != "N/A" else 0.0
                
                # Calculate azimuth directly from GUI values
                azimuth_rad = np.arctan2(gui_east, gui_north)
                az_calc = np.degrees(azimuth_rad)
                if az_calc < 0:
                    az_calc += 360
                    
                print(f"[DEBUG] Using GUI box values for azimuth: atan2({gui_east:.6f}, {gui_north:.6f}) = {az_calc:.2f}°")
                
            except (ValueError, TypeError):
                print(f"[DEBUG] Could not parse GUI box values, falling back to calc_sig")
                az_calc = calculate_signal_azimuth([calc_sig], p_index)  # atan2 result
            back_az_calc = (az_calc + 180) % 360
            p_polarity = dpg.get_value("P Polarity")
            
            # DEBUG: Print polarity and azimuth values
            print(f"[DEBUG] P Polarity read: '{p_polarity}' (type: {type(p_polarity)})")
            print(f"[DEBUG] az_calc: {az_calc:.2f}°, back_az_calc: {back_az_calc:.2f}°")
            
            # atan2 gives direct azimuth, then apply polarity correction:
            # Down motion (compression) = direct azimuth toward epicenter
            # Up motion (dilation) = azimuth + 180° (opposite direction)
            if p_polarity == "down":
                event_az_calc = az_calc  # Compression: direct azimuth
                print(f"[DEBUG] Applied DOWN correction: {event_az_calc:.2f}°")
            elif p_polarity == "up":
                event_az_calc = back_az_calc  # Dilation: opposite direction
                print(f"[DEBUG] Applied UP correction: {event_az_calc:.2f}°")
            else:
                event_az_calc = "-"
                print(f"[DEBUG] No correction applied, polarity unclear: '{p_polarity}'")
                
            dpg.set_value("Azimuth (Calc.)", f"{az_calc:.2f}°")
            dpg.set_value("Back Azimuth (Calc.)", f"{back_az_calc:.2f}°")
            dpg.set_value("Event Azimuth (Calc.)", f"{event_az_calc:.2f}°" if event_az_calc != "-" else "-")
            dpg.split_frame()  # GUI eventlerini işle
        else:
            dpg.set_value("Azimuth (Calc.)", "-")
            dpg.set_value("Back Azimuth (Calc.)", "-")
            dpg.set_value("Event Azimuth (Calc.)", "-")
            dpg.split_frame()
    except Exception as e:
        dpg.set_value("Azimuth (Calc.)", "-")
        dpg.set_value("Back Azimuth (Calc.)", "-")
        dpg.set_value("Event Azimuth (Calc.)", "-")
        dpg.split_frame()

# GUI Layout

dpg.create_context()
dpg.create_viewport(title="Earthquake Signal GUI", width=1553, height=900)

# Register the texture before the main window with a blank map generated by draw_turkey_map
import os
from PIL import Image
MAP_PNG_PATH = os.path.abspath("turkey_map.png")
print(f"[DEBUG] Initial map PNG path: {MAP_PNG_PATH}")
draw_turkey_map(None, None, out_file=MAP_PNG_PATH)
try:
    image = Image.open(MAP_PNG_PATH).convert("RGBA")
    map_data = np.array(image).astype(np.float32) / 255.0
    with dpg.texture_registry():
        dpg.add_dynamic_texture(image.width, image.height, map_data.flatten().tolist(), tag="TurkeyMapTex")
except Exception as e:
    print(f"[DEBUG] Map image load error: {e}")

MAIN_WINDOW_WIDTH = 1550
MAIN_WINDOW_HEIGHT = 950
CHILD_HEIGHT = 950
with dpg.window(label="Earthquake Signal GUI", width=1550, height=950, pos=(0,0)):
    dpg.add_input_text(tag="Selected Folder", default_value=DEFAULT_FOLDER, readonly=True, label="Selected Folder", width=600)
    with dpg.group(horizontal=True):
        # Left child window: File List, plots
        with dpg.child_window(width=600, height=CHILD_HEIGHT):
            dpg.add_text("Files in folder: 0", tag="FileCountText")
            dpg.add_input_text(tag="File Search", hint="Search files...", callback=search_files_callback, width=550, label="")
            dpg.add_listbox(tag="File List", items=[], callback=file_select_callback, width=550, num_items=8, label="")
            for i in range(3):
                # En alttaki Up channel için daha fazla yükseklik
                plot_height = 250 if i == 2 else 200  # Up channel (i=2) için 250, diğerleri için 200
                with dpg.plot(label=f"Channel {i+1} (North)" if i==0 else (f"Channel {i+1} (East)" if i==1 else f"Channel {i+1} (Up)"), height=plot_height, width=550, tag=f"PlotWin##{i}"):
                    dpg.add_plot_axis(dpg.mvXAxis, label="Timestamp", tag=f"PlotAxisX##{i}")
                    with dpg.plot_axis(dpg.mvYAxis, label="Value", tag=f"PlotAxisY##{i}"):
                        dpg.add_line_series([], [], tag=f"Plot##{i}")
                    dpg.add_drag_line(label="P", tag=f"PPickLine##{i}", color=[0,0,255,255], thickness=2, default_value=0.0, callback=sync_cursor_p, show=True)
                    dpg.add_drag_line(label="S", tag=f"SPickLine##{i}", color=[255,0,0,255], thickness=2, default_value=0.0, callback=sync_cursor_s, show=True)
            dpg.add_separator()
        # Middle child window: map, distance, azimuth, and buttons
        with dpg.child_window(width=600, height=CHILD_HEIGHT):
            dpg.add_text("Epicenter (red), Station (blue)")
            dpg.add_image("TurkeyMapTex", width=550, height=350, tag="TurkeyMapImage")
            # Vp/Vs inputları
            with dpg.group(horizontal=True):
                dpg.add_input_float(tag="Vp Value", label="Vp (km/s)", default_value=6.0, min_value=0.1, max_value=20.0, width=100, callback=lambda s,a: vp_vs_changed(), step=0.1)
                dpg.add_input_float(tag="Vs Value", label="Vs (km/s)", default_value=3.0, min_value=0.1, max_value=20.0, width=100, callback=lambda s,a: vp_vs_changed(), step=0.1)
            # Map legend (sağ alt köşe)
            with dpg.group(horizontal=True):
                dpg.add_spacer(width=350)
                with dpg.group(horizontal=True):
                    dpg.add_text("●", color=[255,0,0,255])
                    dpg.add_text(" Epicenter  ")
                    dpg.add_text("●", color=[0,102,255,255])
                    dpg.add_text(" Station")
            dpg.add_separator()
            # Başlıklar
            with dpg.group(horizontal=True):
                dpg.add_text("Original", color=[255,255,255], wrap=150)
                dpg.add_spacer(width=145)
                dpg.add_text("", wrap=120)
                dpg.add_spacer(width=75)
                dpg.add_text("Calculated", color=[255,255,255], wrap=150)
            # Satır bazlı spacer ayarları
            spacer_widths = {
                "Distance": 50,
                "Azimuth": 57,
                "Back Azimuth": 22,
                "Event Azimuth": 15
            }
            for label, org_tag, calc_tag in [
                ("Distance", "Epicenter-Statco Distance", "SP Distance"),
                ("Azimuth", "Azimuth", "Azimuth (Calc.)"),
                ("Back Azimuth", "Back Azimuth", "Back Azimuth (Calc.)"),
                ("Event Azimuth", "Event Azimuth", "Event Azimuth (Calc.)")
            ]:
                with dpg.group(horizontal=True):
                    dpg.add_input_text(tag=org_tag, label="", default_value="-", readonly=True, width=150)
                    dpg.add_spacer(width=20)
                    dpg.add_text(label, color=[255,255,255], wrap=120)
                    dpg.add_spacer(width=spacer_widths[label])
                    dpg.add_input_text(tag=calc_tag, label="", default_value="-", readonly=True, width=150)
                    if label == "Event Azimuth":
                        dpg.add_button(label="Calculate", callback=lambda: file_select_callback(None, None), width=100)
            dpg.add_separator()
            dpg.add_text("P and S Pick Values:")
            # Headers
            with dpg.group(horizontal=True):
                dpg.add_spacer(width=60)  # Space for Pick labels
                dpg.add_text("Time (X)", color=[255,255,255])
                dpg.add_spacer(width=0)
                dpg.add_text("North (Y)", color=[255,100,100])
                dpg.add_spacer(width=0)
                dpg.add_text("East (Y)", color=[100,255,100])
                dpg.add_spacer(width=0)
                dpg.add_text("Up (Y)", color=[100,100,255])
                dpg.add_spacer(width=30)
                dpg.add_text("Action", color=[255,255,255])
            
            # P Pick row
            with dpg.group(horizontal=True):
                dpg.add_text("P Pick:", color=[255,255,100])
                dpg.add_input_text(tag="P Pick", default_value="0.0", readonly=True, width=70, label="")
                dpg.add_input_text(tag="P Y0", default_value="N/A", readonly=True, width=70, label="")
                dpg.add_input_text(tag="P Y1", default_value="N/A", readonly=True, width=70, label="")
                dpg.add_input_text(tag="P Y2", default_value="N/A", readonly=True, width=70, label="")
                dpg.add_button(label="Save Ptime", callback=save_ptime_callback, width=80)
            
            # S Pick row
            with dpg.group(horizontal=True):
                dpg.add_text("S Pick:", color=[255,255,100])
                dpg.add_input_text(tag="S Pick", default_value="0.0", readonly=True, width=70, label="")
                dpg.add_input_text(tag="S Y0", default_value="N/A", readonly=True, width=70, label="")
                dpg.add_input_text(tag="S Y1", default_value="N/A", readonly=True, width=70, label="")
                dpg.add_input_text(tag="S Y2", default_value="N/A", readonly=True, width=70, label="")
                dpg.add_button(label="Save Stime", callback=save_stime_callback, width=80)
            dpg.add_separator()
            dpg.add_input_text(tag="P Polarity", label="P Polarity", default_value="-", readonly=True, width=120)
            dpg.add_separator()
            dpg.add_text("Export all files to Excel:")
            dpg.add_button(label="Save Excel", callback=save_excel_callback, width=200)
            dpg.add_progress_bar(tag="Excel Progress", default_value=0.0, width=200)
            dpg.add_text("", tag="Progress Text")
            dpg.add_separator()
        # Right child window: attributes
        with dpg.child_window(width=350, height=CHILD_HEIGHT):
            dpg.add_text("Available attributes in mat['EQ'].anEQ:")
            dpg.add_input_text(tag="anEQ Attributes List", multiline=True, readonly=True, width=310, height=CHILD_HEIGHT-50, default_value="-")

with dpg.handler_registry():
    pass  # no render callback

# Vp/Vs değişince dosyayı ve hesaplamaları güncelleyen fonksiyon
def vp_vs_changed():
    file_select_callback(None, None)

# Initialize file list after GUI is fully ready
initialize_file_list()

dpg.setup_dearpygui()
dpg.show_viewport()

dpg.start_dearpygui()
dpg.destroy_context() 