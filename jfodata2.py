# -*- coding: utf-8 -*-
"""
Utility functions for data plotting and processing.

Originally developed by Jaime for routine analysis in jfodata package.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate

# Optional: customize global matplotlib settings
plt.rcParams.update({
    'axes.linewidth': 1.4,
    'figure.figsize': (12, 9),
    'font.size': 32,
    'legend.fontsize': 'large',
    'figure.titlesize': 'medium',
    'axes.labelsize': 32,
    'savefig.bbox': 'tight',
})


def interpolate_dataframe(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    x_min: float,
    x_max: float,
    n_points: int = 1000,
    smoothing: float = 1.0,
    spline_degree: int = 3
):
    """
    Interpolates y values using a univariate spline over a specified x range.

    Parameters:
        df (pd.DataFrame): Input dataframe.
        x_col (str): Name of x column.
        y_col (str): Name of y column.
        x_min (float): Minimum x value for interpolation.
        x_max (float): Maximum x value for interpolation.
        n_points (int): Number of interpolation points.
        smoothing (float): Smoothing factor.
        spline_degree (int): Degree of the spline (default 3 = cubic).

    Returns:
        tuple: (spline function, interpolated x values)
    """
    df = df.copy()
    df = df.dropna(axis=1, how='all').drop_duplicates(subset=x_col)
    df = df.sort_values(by=x_col)
    df = df[(df[x_col] > x_min) & (df[x_col] < x_max)]

    x_interp = np.linspace(x_min, x_max, n_points)
    spline = interpolate.UnivariateSpline(df[x_col], df[y_col],
                                          k=spline_degree,
                                          s=smoothing,
                                          ext=3)
    return spline, x_interp


def list_files_by_extension(folder_path: str, extension: str, show: bool = True):
    """
    Lists all files with a given extension in the specified folder.

    Parameters:
        folder_path (str): Path to the folder.
        extension (str): File extension (e.g., 'txt', 'csv').
        show (bool): Whether to print the list.

    Returns:
        tuple: (list of filenames, count)
    """
    os.chdir(folder_path)
    files = [f for f in os.listdir(".") if f.endswith("." + extension)]
    if show:
        for i, file in enumerate(files):
            print(f"{i} {file}")
    return files, len(files)


def import_data_from_folder(folder_path: str, extension: str, index: int,
                            separator: str = '\t', skip_rows: int = 1):
    """
    Imports data from a selected file in a folder.

    Parameters:
        folder_path (str): Path to the folder.
        extension (str): File extension.
        index (int): Index of the file to import.
        separator (str): Column separator.
        skip_rows (int): Rows to skip at the start.

    Returns:
        pd.DataFrame: Loaded data.
    """
    os.chdir(folder_path)
    files = [f for f in os.listdir(".") if f.endswith("." + extension)]
    return pd.read_csv(files[index], sep=separator, skiprows=skip_rows)


def create_resistance_vs_temperature(df, temp_col: str, resistance_col: str):
    """
    Creates a DataFrame with calculated resistivity and conductivity from temperature and resistance.

    Returns:
        pd.DataFrame
    """
    result = pd.DataFrame()
    result['T'] = df[temp_col]
    result['R'] = df[resistance_col]
    result['rho'] = (0.05 / 1.221) * df[resistance_col]
    result['sigma'] = 1 / result['rho']
    return result


def create_magnetoresistance(df, temp_col: str, resistance_col: str,
                             field_col: str, field_rate: float):
    """
    Processes magnetoresistance data with calculated field and conductivity.

    Returns:
        pd.DataFrame
    """
    result = create_resistance_vs_temperature(df, temp_col, resistance_col)
    result['Field'] = field_rate * df[field_col]
    return result


def create_mr_from_voltage(df, temp_col: str, voltage_col: str,
                           field_col: str, current_mA: float):
    """
    Processes magnetoresistance data from voltage measurements.

    Returns:
        pd.DataFrame
    """
    result = pd.DataFrame()
    result['T (K)'] = df[temp_col]
    result['R (Ohm)'] = df[voltage_col] / (current_mA / 1000)
    result['Field (T)'] = df[field_col] / 10000
    return result
    
 def find_closest_index(df: pd.DataFrame, column: str, value: float) -> int:
    """
    Returns the index of the row in the DataFrame where the value in the specified column is closest to the given value.
    
    Parameters:
        df (pd.DataFrame): DataFrame to search.
        column (str): Column name to check.
        value (float): Target value to find closest match.
        
    Returns:
        int: Index of the closest row.
    """
    return (np.abs(df[column] - value)).idxmin()

def separate_qv_curves(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extracts and processes two dataframes from the input dataframe from quantum desing files, corresponding to
    two measurement channels (Ch1 and Ch2). Converts units, renames columns for
    standardization, drops rows with missing data, and resets indices.
    
    Parameters:
        df (pd.DataFrame): Raw input dataframe with original measurement columns.
        
    Returns:
        tuple: Two pd.DataFrames (df1, df2) corresponding to channel 1 and channel 2.
               If extraction fails for a channel, returns empty DataFrame for that channel.
    """
    # Helper to safely extract channel data
    def extract_channel_data(channel_num: int) -> pd.DataFrame:
        try:
            ch = str(channel_num)
            data = {
                'T(K)': df['Temperature (K)'],
                'B(T)': round(df['Field (Oe)'] / 10000, 4 if channel_num == 1 else 2),
                'Angle (deg)': df['Sample Position (deg)'],
                'V(V)': df[f'In Phase Voltage Ampl Ch{ch} (V)'],
                '4V(V)': df[f'Quadrature Voltage Ch{ch} (V)'],
                'R(ohms)': df[f'Resistance Ch{ch} (Ohms)'],
                'Phase(deg)': df[f'Phase Angle Ch{ch} (deg)'],
                'I(A)': df[f'AC Current Ch{ch} (mA)'] / 1000,
                'f(Hz)': df[f'Frequency Ch{ch} (Hz)'],
                '2 Harm': df[f'2nd Harmonic Ch{ch} (dB)'],
                'W(rad/s)': 2 * np.pi * df[f'Frequency Ch{ch} (Hz)'],
            }
            if channel_num == 1:
                # Additional error column only for channel 1
                data['R1_erro(ohms)'] = df['Resistance Std. Dev. Ch1 (Ohms)']
            
            df_channel = pd.DataFrame(data)
            df_channel.dropna(inplace=True)
            df_channel.reset_index(drop=True, inplace=True)
            return df_channel
        
        except KeyError:
            # If any column is missing, return empty DataFrame
            return pd.DataFrame()

    df1 = extract_channel_data(1)
    df2 = extract_channel_data(2)

    return df1, df2
    
    
def detect_data_type(df):
    """
    Determine the type of measurement data in the DataFrame: 
    Magnetoresistance (MR), Angular scan, or Resistance vs Temperature (RvsT).

    Parameters:
        df (pd.DataFrame): Input dataframe with measurement data.

    Returns:
        tuple: (mode, ch1_present, ch2_present)
            mode (int):
                0 -> RvsT (Resistance vs Temperature)
                1 -> MR (Magnetoresistance)
                2 -> Angular scan
                3 -> Unknown / Error
            ch1_present (int): 1 if Resistance Ch1 data is present, else 0
            ch2_present (int): 1 if Resistance Ch2 data is present, else 0
    """
    ch1_present = 0
    ch2_present = 0

    try:
        # Check variance of Temperature and Field columns to detect mode
        if round(df['Temperature (K)'].var()) != 0:
            mode = 0  # RvsT
        elif round(df['Field (Oe)'].var()) != 0:
            mode = 1  # MR
        else:
            mode = 2  # Angular scan

        # Check if Resistance columns have any non-NaN values
        if df['Resistance Ch1 (Ohms)'].dropna().size > 0:
            ch1_present = 1
        if df['Resistance Ch2 (Ohms)'].dropna().size > 0:
            ch2_present = 1

    except (KeyError, AttributeError):
        # If expected columns are missing or df is invalid
        mode = 3

    return mode, ch1_present, ch2_present    
