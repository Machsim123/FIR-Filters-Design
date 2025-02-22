# **FIR Filters Design**  

## 📌 **Project Overview**  
This project focuses on designing **Finite Impulse Response (FIR) filters** using **optimization techniques** in MATLAB. It explores two primary methods:

- **Least Squares Method (L2 norm - MMSE)**: Minimizes the mean squared error between the ideal and actual filter response.  
- **Chebyshev Method (L∞ norm - Minimax)**: Uses the **Parks-McClellan algorithm** to achieve an equiripple response.  

The project also evaluates the **frequency response**, **phase characteristics**, and **stopband attenuation** of different filters with various orders.  

---  

## 📂 **Features & Functionality**  
✅ **FIR Filter Design:** Implements optimal FIR filters using `firls` (Least Squares) and `firpm`/`remez` (Chebyshev).  
✅ **Frequency Response Analysis:** Uses `freqz` to visualize magnitude and phase response.  
✅ **Comparison of Filter Designs:** Analyzes differences in **error norms, ripple behavior, and transition width**.  
✅ **Weighting Function Implementation:** Enhances stopband attenuation using weighted optimization.  
✅ **MATLAB Scripting:** Fully automated design process with adjustable filter order (M = 10 to 100).  

---  

## 📈 **Results & Analysis**  
- **Least Squares Method** results in a **smoother response**, with a **small mean squared error** but wider transition bands.  
- **Chebyshev Method** achieves **better stopband attenuation** but introduces **equiripple behavior** in the passband.  
- Weighting functions significantly **enhance stopband attenuation** by penalizing errors in certain frequency regions.  

---  

## 🚀 **Future Improvements**  
🔹 Extend to **IIR filter design** using Butterworth, Chebyshev, and Elliptic filters.  
🔹 Implement **adaptive filtering** for real-time signal processing applications.  
🔹 Convert MATLAB implementation to **Python (SciPy, NumPy) for wider accessibility**.  

---  

## 🏆 **Credits**  
Developed by **Maxim Vaculenco**  
📧 Contact: maxim.vaculenco.jr@gmail.com  
