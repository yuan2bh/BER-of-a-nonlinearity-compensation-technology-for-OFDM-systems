# BER-of-a-nonlinearity-compensation-technology-for-OFDM-systems
A nonlinear Hammerstein system parameter separation technique was used for nonlinear and linear coefficient separation, which is not perfect.   

The abstract of the relevant published paper:   
Aiming at the problem of communication performance of OFDM systems when the Power Amplifier (PA) nonlinearity is strong, and based on a kind of combined estimation technology for nonlinearity of transmitters and wireless channels, a nonlinear channel estimation and equalization technique of OFDM system is proposed based on training sequences. Firstly, the transmitter nonlinearity and unit impulse response of wireless channel are estimated with least squares (LS) algorithm, then the compensation of wireless channel and transmitter nonlinearity are done successively. Simulation results show that the proposed method may approach the result of OFDM wireless communication systems without nonlinearity under the perfect compensation of wireless channe1.     
Key words: OFDM; nonlinear channel; channel estimation; channel compensation; LS   

The published paper:  
袁红林, 刘成云, 包志华. 一种正交频分复用系统的非线性补偿技术[J]. 电视技术, 2016, 40(1): 87-90, 94.     
YUAN H.L, LIU C. Y, BAO Z. H. Nonlinearity compensation technology for OFDM systems[J]. Video engineering, 2016, 40(1): 87-90. (in Chinese)     

The entry is sepaComp_main.m which can run successfully in MATLAB R2013b.   
One of the BER comparisons between our proposed method and the theoretical value of Rayleigh fading channel is as follows:    
![image](https://github.com/user-attachments/assets/3bb59765-7653-47d3-9eda-084145ebe9c9) 

The BER of our proposed method is worse than that of the theoretical Rayleigh fading channel for a number of reasons, including the fact that the theoretical values do not take into account nonlinear effects, as well as a number of other research-worthy issues.  
