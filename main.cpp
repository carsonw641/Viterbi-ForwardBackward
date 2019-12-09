#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stack>
#include <string>
#include <vector>

double logVal(double val) {
    if (val > 0) {
        return std::log10(val);
    }

    return std::numeric_limits<double>::lowest();
}

void prettyPrint(std::vector<std::vector<int>>& traceBack, double endingValue, int endState){
    std::stack<int> hist;

    int iter = endState;
    int x = traceBack[endState].size()-1;

    while(iter != -1){
        hist.push(traceBack[iter][x]);
        iter = hist.top();
        x--;
    }

    while (hist.size() > 0){
        int state = hist.top();

        if (state < 0){
            std::cout << "Start -> ";
        }else {
            std::cout << state << " -> ";
        }

        hist.pop();
    }

    std::cout << endState << " ending value : " << endingValue << std::endl;

}

void viterbi(std::vector<std::vector<double>>& sMatrix,
            std::vector<std::vector<double>>& eMatrix,
            std::vector<double>& initalStates,
            std::vector<double>& endingStates,
            std::vector<std::string>& sequence,
            std::map<std::string, int>& titles) {

    std::vector<double> output (sMatrix.size());
    std::vector<std::vector<int>> traceBack (sMatrix.size(), std::vector<int>(sequence.size(), -1));

    // pre-calculate the inital emission with initial state probabilities 
    for (int i = 0; i < initalStates.size(); i++){
        output[i] = logVal(initalStates[i] * eMatrix[i][titles[sequence[0]]]);
    }
    
    for (int seq = 1; seq < sequence.size(); seq++){
        std::vector<double> tempOutput (sMatrix.size());

        // the next 
        for (int col = 0; col < sMatrix.size(); col++){
            double maxVal = std::numeric_limits<double>::lowest();
            int prevState = -1;
            
            double emitProb = eMatrix[col][titles[sequence[seq]]];

            // determining which state to transfer from
            if (emitProb > 0) {
                for (int row = 0; row < sMatrix.size(); row++){
                    if (output[row] > std::numeric_limits<double>::lowest() && sMatrix[row][col] > 0){
                        double val = output[row] + logVal(emitProb * sMatrix[row][col]);

                        if (val > maxVal){
                            maxVal = val;
                            prevState = row;
                        }
                    }
                }
            }

            traceBack[col][seq] = prevState;
            tempOutput[col]  = maxVal;
        }

        output = tempOutput;
    }

    for (int i = 0; i < endingStates.size(); i++) {
        output[i] = output[i] + logVal(endingStates[i]);
        
        prettyPrint(traceBack, output[i], i);
    }

}

double forward(std::vector<std::vector<double>>& sMatrix,
            std::vector<std::vector<double>>& eMatrix,
            std::vector<double>& initalStates,
            std::vector<double>& endingStates,
            std::vector<std::string>& sequence,
            std::map<std::string, int>& titles,
            int e, int s){

    std::vector<double> output (sMatrix.size());

    // pre-calculate the inital emission with initial state probabilities 
    for (int i = 0; i < initalStates.size(); i++){
        output[i] = initalStates[i] * eMatrix[i][titles[sequence[0]]];
    }


    for (int seq = 1; seq < e; seq++){
        if (seq == sequence.size()){
            double total = 0;
            for (int i = 0; i < output.size(); i++){
                total += output[i]*endingStates[i];
            }

            return total;
        }
        
        std::vector<double> tempOutput (sMatrix.size());

        // the next 
        for (int col = 0; col < sMatrix.size(); col++){
            double sum = 0;
            
            double emitProb = eMatrix[col][titles[sequence[seq]]];

            // determining which state to transfer from
            if (emitProb > 0) {
                for (int row = 0; row < sMatrix.size(); row++){
                    if (output[row] > std::numeric_limits<double>::lowest() && sMatrix[row][col] > 0){
                        sum += output[row]*emitProb*sMatrix[row][col];
                    }
                }
            }

            if (sum != 0) {
                tempOutput[col] = sum;
            }else {
                tempOutput[col] = std::numeric_limits<double>::lowest();
            }
            
        }

        output = tempOutput;
    }

    return output[s-1];
}

double backward(std::vector<std::vector<double>>& sMatrix,
            std::vector<std::vector<double>>& eMatrix,
            std::vector<double>& initalStates,
            std::vector<double>& endingStates,
            std::vector<std::string>& sequence,
            std::map<std::string, int>& titles,
            int e, int s){
    
    std::vector<double> output (sMatrix.size());

    for (int i = 0; i < endingStates.size(); i++){
        output[i] = endingStates[i];
    }

    for (int seq = sequence.size()-1; seq >= e; seq--){
        std::vector<double> tempOutput (sMatrix.size());

        //state transitioning from
        for (int col = 0; col < sMatrix.size(); col++){
            double sum = 0;
            // state transitioning too
            for (int row = 0; row < sMatrix.size(); row++){
                if (output[row] != std::numeric_limits<double>::lowest()){
                    double transitionVal = sMatrix[col][row]*eMatrix[row][titles[sequence[seq]]]*output[row];
                    sum += transitionVal;
                }
            }
            
            if (sum != 0) {
                tempOutput[col] = sum;
            }else {
                tempOutput[col] = std::numeric_limits<double>::lowest();
            }
            
        }

        output = tempOutput;
    }

    return output[s-1];
}

int main () {

    std::ifstream fileInput("config", std::ios::binary);

    if (!fileInput.is_open()){
        std::cout << "Failed to open config file \n";
    }else {

        /* BEGIN INPUT */

        int states, emissions;
        std::string input;

        // Read input for states and emissions
        fileInput >> states >> emissions;

        std::vector<std::vector<double>> sMatrix (states, std::vector<double>(states));
        std::vector<std::vector<double>> eMatrix (states, std::vector<double>(emissions));

        std::vector<double> initalStates (states);
        std::vector<double> endingStates (states);

        std::map<std::string, int> emissionTitles;
        std::vector<std::string> sequence;

        // Read initial state probabilities
        for (int i = 0; i < states; i++){
            fileInput >> initalStates[i];
        }

        // Read ending state probabilities
        for (int i = 0; i < states; i++){
            fileInput >> endingStates[i];
        }

        // Read state transition matrix
        for (int i = 0; i < states; i++){
            for (int j = 0; j < states; j++){
                fileInput >> sMatrix[i][j];
            }
        }

        // Read emission emissionTitles
        fileInput >> input;
        input = input.substr(1, input.size()-2);

        size_t pos = 0;
        int iter = 0;
        while ((pos = input.find(',')) != std::string::npos) {
            emissionTitles[input.substr(0, pos)] = iter;
            input.erase(0, pos+1);
            iter++;
        }
        emissionTitles[input.substr(0, pos)] = iter;

        // Read emission values
        for (int i = 0; i < states; i++){
            for (int j = 0; j < emissions; j++){
                fileInput >> eMatrix[i][j];
            }
        }

        // Read sequence
        fileInput >> input;
        input = input.substr(1, input.size()-2);

        pos = 0;
        while ((pos = input.find(',')) != std::string::npos) {
            sequence.push_back(input.substr(0, pos));
            input.erase(0, pos+1);
        }
        sequence.push_back(input.substr(0, pos));

        //  Read optional query
        int e = -1, s = -1;

        if (!fileInput.eof()){
            fileInput >> input;
            input = input.substr(1, input.size()-2);
            e = std::stoi(input.substr(0, 1));
            s = std::stoi(input.substr(2, 1));
        }

        /* END INPUT */
        double fwd, bkw;
        viterbi(sMatrix, eMatrix, initalStates, endingStates, sequence, emissionTitles);
        std::cout << "Forward: " << (fwd = forward(sMatrix, eMatrix, initalStates, endingStates, sequence, emissionTitles, e, s)) << std::endl;
        std::cout << "Backward: " << (bkw = backward(sMatrix, eMatrix, initalStates, endingStates, sequence, emissionTitles, e, s)) << std::endl;
        std::cout << "P(x): " << ((fwd*bkw)/forward(sMatrix, eMatrix, initalStates, endingStates, sequence, emissionTitles, sequence.size()+1, 0)) << std::endl;

    }

    fileInput.close();

    return 0;
}