# OFDM_Audio_Transmission_Project

Return to [Joy's website](https://joyyeh2002.github.io/engineering.html)

**EE-442 Fall 2022: Wireless Device Algorithms**

**Final Project: Implementing wireless transmission of binary data through an audio system and MATLAB OFDM transmitter and receiver.**

# Introduction
Wideband single carrier systems require very complicated equalization and channel estimation schemes due to frequency selective fading. Orthogonal frequency division multiplexing (OFDM) is a way to convert a wideband channel into many intersymbol interference free narrowband channels, for which equalization and channel estimation are straightforward.


![](/info/UCI_Repository.PNG)

This task is one of professor Hung-Yi Lee’s undergraduate Machine Learning course at National Taiwan University. Assignment instructions can be found [here](https://drive.google.com/file/d/1DOqby_K_9TK2QVUr-C2qYg0amZkRu21U/view). The dataset is revised so make the features easier to manipulate. Relevant data is located in the ‘data’ folder of this GigHub repository.

The dataset’s input, **x_train**, has vectors with 510 dimensions. The training outputs, **y_train**, are binary values of 0 (income below $50k) or 1 (income above $50k). There are 48842 instances of labeled data. Because we have ample amounts of data and simple outputs, we use supervised learning to train a model with logistic regression. Simple implementations of **Adagrad** and **mini-batch training** are also made.

# Inputs and Outputs
The inputs and adjustable parameters are:
- Epoch size
- Initial learning step for Adagrad
- Percentage split for the training and validation sets

The outputs of the model are:

- A graph of the loss values of the training set and validation sets over time
- A graph of the accuracy % of the model over time
- A .csv file that predicts the y values of the testing data

# Workflow
1. Import and clean up the data files: x_train, y_train, and x_test
2. Split up the training and validation sets
3. Normalize the training sets so that for all dimensions of x,  means = 0 and variances = 1
4. Compute loss values of each set of parameters passed in
5. Write functions for shuffling, sigmoid(), creating a proposed function _f(x, w, b), prediction, accuracy, cross entropy, and gradient descent
6. Implement mini-batch and Adagrad
7. Print and update the parameters in the console
8. Plot the graphs mentioned above
9. Generate prediction results and write to the output file, y_test.
![](/info/Code_Running.gif)


# The Statistical Model
This is a logistic regression model that assumes the dimensions to follow Gaussian Distributions. One co-variance matrix is shared between all input dimensions so that there's a linear boundary for the binary classification (whether the output probability is > 0.5). 


# Results
![](/outputs/loss_Gradient.png)

![](/outputs/Accuracy_Gradient.png)

Return to [Joy's website](https://joyyeh2002.github.io/compsci.html)
