""" (Se her Steinar)
Krever python3. Skriv deretter i en terminal:
pip install numpy
pip install matplotlib
Da er du klar! Kan være lurt å ikke kjøre via Visual Studio, men bruke en annen IDE som gjør mindre hokus pokus

Bruk: Fyll ut "filename". Du har også to flag på plotTemperatureEvolution for hvordan du vil dataen skal plottes
"""
#imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import json
from matplotlib.lines import Line2D

# --------------------------- INPUT FILE ---------------------------
# Skriv full fil path ti len .Temp fil, feks har jeg på min maskin:
filename = 'C:/Users/fredr/Oliasoft/welltemp/test/Fredrik_tests/simple_design_prooduction-C1Production16casing.Temp'

#region --------------------------- READ FILE & TEMPS ---------------------------
def readJSONFile(filename):
    with open(filename, 'r') as file:
        json_data = json.load(file)
    return json_data

def fetchTemps(jsonTemp) -> list[np.ndarray]:
    #Get all temps and times
    timeKeys = list(jsonTemp["Time"])
    temperatures = [jsonTemp["Time"][timeKey]["Temperature"] for timeKey in timeKeys]
    return [np.array(timeKeys).astype(float), np.array(temperatures)] #Returns times and the heat map at each time

#endregion Read file -----------------------------------------------------

#region --------------------------- PLOTTING FUNCTION  ---------------------------
def _pltgrid(ax, xBounds, yBounds, heatMap, vmin, vmax):
    # Update the heatmap data only, without clearing the axis
    cax = ax.pcolormesh(xBounds, yBounds, heatMap, shading='auto', cmap='turbo', vmin = vmin, vmax = vmax) #viridis
    return cax

def format_time(seconds) -> str:
    """Convert seconds to a formatted string like '2D 4H 45M'."""
    days = seconds // (24 * 3600)
    hours = (seconds % (24 * 3600)) // 3600
    minutes = (seconds % 3600) // 60
    return f"{int(days)}D {int(hours)}H {int(minutes)}M"

def format_time_seconds(seconds) -> str:
    return f"{(seconds)}S"

def plotTemperatureEvolution(xBounds: np.ndarray, yBounds: np.ndarray, heatEvolution: np.ndarray, showgrid=True, times = [], useGlobalTempSpan = True, label = "Temperature"):
    """Plots the heat map on an uneven grid specified by xBounds & yBounds with a time slider. 
    Args:
        xBounds (np.ndarray): NX + 1 elements. (monotonically increasing / decreasing)
        yBounds (np.ndarray): NY + 1 elements. (monotonically increasing / decreasing)
        heatEvolution (np.ndarray): Shape (NT, NY, NX). First index is time, second is Y, third is X
        showgrid (bool, optional): Show the grid lines. Defaults to True.
        useGlobalTempSpan (bool, optional): If enabled, the coloring will be done based on max/min of the data for the entire heat evolution.
            It gives a more true feeling of what is hot and cold since it won't always use the full spectrum. Downside is that it can make  
            it hard to read nuances if the temp san is much greater at some time series' than others. Defaults to True.
    """
    # Set up figure and initial plot
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.invert_yaxis()
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)  # Make room for the slider

    # Initial plot
	#Set span of data based on total min an max
    if useGlobalTempSpan:
        global_min = np.min(heatEvolution); global_max = np.max(heatEvolution)
    else: global_min = global_max = None
    cax = _pltgrid(ax, xBounds, yBounds, heatEvolution[0], global_min, global_max)
    cbar = fig.colorbar(cax, ax=ax, label=label)

    # Draw grid lines as Line2D objects, only if showgrid is True
    if showgrid:
        # Create vertical grid lines
        for x in xBounds:
            line = Line2D([x, x], [yBounds[0], yBounds[-1]], color='black', linestyle='--', linewidth=0.5)
            ax.add_line(line)
        
        # Create horizontal grid lines
        for y in yBounds:
            line = Line2D([xBounds[0], xBounds[-1]], [y, y], color='black', linestyle='--', linewidth=0.5)
            ax.add_line(line)

    # Slider for selecting time
    ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03], facecolor="lightgray")  # type: ignore # Position and size of slider
    slider = Slider(ax_slider, "Time", 0, len(heatEvolution) - 1, valinit=0, valfmt="%0.0f")
    # Update the slider label with formatted time
    if (len(times) > 0): slider.valtext.set_text(format_time(times[0])) # type: ignore

    # Add a text box to display temperature with a dummy initial value
    text = ax.text(0.02, 0.95, label + ": --.- at X=--.-, Y=--.-", 
                   transform=ax.transAxes, fontsize=12, verticalalignment='top')

    # Function to update the heatmap when the slider is changed
    def update(val):
        idx = int(slider.val)
        # Update the heatmap without clearing the grid lines
        cax = _pltgrid(ax, xBounds, yBounds, heatEvolution[idx], global_min, global_max)
        # Update the slider label with the correct time format
        if (len(times) > 0): slider.valtext.set_text(format_time(times[idx])) # type: ignore
        # Update colorbar if you don't use global temp span
        if (not useGlobalTempSpan): cbar.update_normal(cax)
        fig.canvas.draw_idle()

    # Function to update the displayed temperature on hover
    def on_hover(event):
        if event.inaxes == ax:
            x, y = event.xdata, event.ydata
            # Find the indices of the corresponding grid cell
            if x is not None and y is not None:
                col = np.searchsorted(xBounds, x) - 1
                row = np.searchsorted(yBounds, y) - 1
                if 0 <= col < heatEvolution.shape[2] and 0 <= row < heatEvolution.shape[1]:
                    idx = int(slider.val)
                    temp = heatEvolution[idx, row, col]
                    text.set_text(label + f": {temp:.2f} at X={x:.2f}, Y={y:.2f}")
                else:
                    text.set_text(label + f": --.- at X={x:.2f}, Y={y:.2f}.")
            fig.canvas.draw_idle()

    # Connect the hover event to the on_hover function
    fig.canvas.mpl_connect("motion_notify_event", on_hover)
    # Connect the update function to the slider
    slider.on_changed(update)
    # Show the plot
    plt.show()

#endregion Plotting function -----------------------------------------------------

# Read and parse .. 
jsonTemp = readJSONFile(filename)
print('\nFile loaded with Simulation ID: ', jsonTemp['SimID'])
[times, heatEvolution] = fetchTemps(jsonTemp)
#Prepare grid
zBounds = np.array(jsonTemp["zBound"])
rBounds = np.array(jsonTemp["rBound"])

# Heat evolution plot (Se her Steinar)
#Noen ting du kan velge: showgrid = True gjør at man ser gridet som ble brukt i plottet
#useGlobalTempSpan gjør at heatmapet bruker farger ut ifra maks og min temperatur for hele dataserien. Hvis du skrur den av så vil
#	den bruke hele temperaturskalaen for hvert stillbilde i stedet. Det gjør det lettere å se forskjeller i temperatur for en gitt tid, men
#	til gjengjelg kan det være misvisende på hva som virkelig er kaldt og virkelig er varmt.
plotTemperatureEvolution(rBounds, zBounds, heatEvolution, showgrid=True, times = times, useGlobalTempSpan=True)