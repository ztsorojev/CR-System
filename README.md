# Cognitive Radio System

The goal of this project was to develop a **MATLAB code** to determine the modulation type of a received signal under AWGN channel in a cognitive radio (CR) system. The signal can be multi-carrier or single carrier and from any of these 5 types of modulations: OFDM, 4-QAM, 16-QAM, 2-PAM or 4-PAM. 

In order to achieve our goal, we implement a classification tree. This tree classifies the signal step by step. We start with the most general classification: multi-carrier vs single carrier. Then, in the case of a single carrier, the next step is to determine the modulation type: M-QAM or M-PAM. This part is done by comparing the feature vector of our signal with the ideal features of QAM and PAM. Once we know that, we extract the symbols from the down-converted version of our received signal and we use their distribution to determine the modulation level. This report will explain in more details how this tree is implemented.


## Files

* The code for the entire classifier is in the file "full_classifier.m".
* The code for each individual block is in the .m file with its name. 
* *Project_report.pdf* contains mathematical explanations about the model, as well as discussions about my results.
