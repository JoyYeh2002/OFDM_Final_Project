# OFDM Audio Transmission Project

Return to [Joy's website](https://joyyeh2002.github.io/engineering.html)

**EE-442 Fall 2022: Wireless Device Algorithms**

**Final Project: Implementing wireless transmission of binary data through an audio system and MATLAB OFDM transmitter and receiver.**

# Introduction
Wideband single carrier systems require very complicated equalization and channel estimation schemes due to frequency selective fading. Orthogonal frequency division multiplexing (OFDM) is a way to convert a wideband channel into many intersymbol interference free narrowband channels, for which equalization and channel estimation are straightforward.

![](/info/OFDM_comparisons.PNG)

This is the final project of a graduate-level wireless communications course at EPFL, Switzerland. Assignment instructions can be found [here](OFDM_Project_Instructions.pdf). The main task is to use **orthogonal frequency division multiplexing** to organize input binary data packets, transmit the data blocks with **block type** and **comb type** transmissions, then use cyclic prefixs and **phase correction** to recover the input data.

![](/info/block_diagram.PNG)

The channels are built from a physical audio system with speakers and a microphone. Therefore, the system is subject to noise and fading effects. With OFDM, we are able to achieve high recovery rates under low SNR (signal to noise ratio). In the MATLAB simulations, we also conducted significant amount of [testing](https://github.com/JoyYeh2002/OFDM_Final_Project/tree/main/Chu_Miao_Yeh_OFDM_Project_Code) to fully evaluate our implementation.

# Main Result Highlights

![](/info/demo01.PNG)
![](/info/demo02.PNG)
![](/info/demo03.PNG)
![](/info/demo04.PNG)
![](/info/demo05.PNG)
![](/info/demo06.PNG)
![](/info/demo07.PNG)
![](/info/demo08.PNG)
![](/info/demo09.PNG)

- The [full demonstration slides](https://github.com/JoyYeh2002/OFDM_Final_Project/blob/main/Chu_Miao_Yeh_OFDM_Project_Presentation_PDF.pdf)
- The [full project report](https://github.com/JoyYeh2002/OFDM_Final_Project/blob/main/Chu_Miao_Yeh_OFDM_Project_Report.pdf)

Return to [Joy's website](https://joyyeh2002.github.io/engineering.html) 
