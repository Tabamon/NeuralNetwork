//Lorenzo Toscano

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>

struct Connection
{
    double weight;
    double deltaWeight;
};

class Neuron;
typedef std::vector<Neuron> Layer;

class Neuron{
    public:
        Neuron(unsigned int numOfOutputs, unsigned int myIndex);
        void setOutputVal (double value);
        void neuronFeedForward(const Layer &prevLayer);
        double getOutputVal (void) const;
        void calcOutputGradient(double targetVal);
        void calcHiddenGradient(const Layer &nextLayer);
        void updateInputWeights(Layer &prevLayer);

    private:
        static double eta;
        static double alpha;
        static double randomWeight(void) {return rand() / double(RAND_MAX);}
        double neuronOutput;
        std::vector<Connection> outputWeights;
        unsigned int myIndex;
        static double transformationFunction(double x);
        static double transformationFunctionDerivative(double x);
        double gradient;
        double sumDOW(const Layer & nextLayer) const;
};

double Neuron::alpha = 0.5;
double Neuron::eta = 0.15;

void Neuron::updateInputWeights(Layer &prevLayer){
    for(unsigned int neuronCounter =0; neuronCounter<prevLayer.size(); neuronCounter++){
        Neuron &neuron = prevLayer[neuronCounter];
        double oldDeltaWeight = neuron.outputWeights[myIndex].deltaWeight;
        double newDeltaWeight = eta * neuron.getOutputVal() * gradient + alpha * oldDeltaWeight;
        neuron.outputWeights[myIndex].deltaWeight = newDeltaWeight;
        neuron.outputWeights[myIndex].weight += newDeltaWeight;
    }
}

double Neuron::sumDOW(const Layer & nextLayer) const{
    double sum = 0.0;
    for(unsigned int neuronCounter =0; neuronCounter<nextLayer.size(); neuronCounter++){
        sum += outputWeights[neuronCounter].weight * nextLayer[neuronCounter].gradient;
    }
    return sum;
}

void Neuron::calcHiddenGradient(const Layer &nextLayer){
    double dow = sumDOW(nextLayer);
    gradient = dow * transformationFunctionDerivative(neuronOutput);
}

void Neuron::calcOutputGradient(double targetVal){
    double delta = targetVal - neuronOutput;
    gradient = delta * transformationFunctionDerivative(neuronOutput);
}

double Neuron::transformationFunction(double x){
    return tanh(x);
}

double Neuron::transformationFunctionDerivative(double x){
    return 1-x*x;
}

Neuron::Neuron(unsigned int numOfOutputs, unsigned int myIndex){
    for (unsigned c = 0; c < numOfOutputs; c++){
        outputWeights.push_back(Connection());
        outputWeights.back().weight = randomWeight();
    }
    this -> myIndex = myIndex;
}

void Neuron::setOutputVal(double value){
    neuronOutput = value;
}

double Neuron::getOutputVal(void) const{
    return neuronOutput;
}

void Neuron::neuronFeedForward(const Layer &prevLayer){
    double sum = 0.0;

    for (unsigned int neuronCounter = 0; neuronCounter < prevLayer.size(); neuronCounter++){
        sum += prevLayer[neuronCounter].getOutputVal()*prevLayer[neuronCounter].outputWeights[myIndex].weight;
    }
    neuronOutput = transformationFunction(sum);
}

//Class Net----------------------------------------

class Net{
    public:
        Net(const std::vector<unsigned int> &topology);
        void feedForward(const std::vector<double> &inputVals);
        void backProp(const std::vector<double> &targetVals);
        void getResults(std::vector<double> &resultVals) const;
   
    private:
        std::vector<Layer> layers; //matrix of neurons
        double error;
        double recentAverageError;
        double recentAverageSmoothigFactor;
};

void Net::getResults(std::vector<double> &resultVals) const{
    resultVals.clear();

    for (unsigned int neuronCounter = 0; neuronCounter < layers.back().size() - 1; neuronCounter++){
        resultVals.push_back(layers.back()[neuronCounter].getOutputVal());
    }
}

Net::Net(const std::vector<unsigned> &topology){

    unsigned int numOfLayers = topology.size();
    for (unsigned int layerCounter = 0; layerCounter<numOfLayers; ++layerCounter){
        layers.push_back(Layer());
        //unsigned numOfOutputs = layerCounter == topology.size() - 1 ? 0 : topology[layerCounter + 1];
        unsigned int numOfOutputs;
        if (layerCounter == topology.size() - 1) numOfOutputs = 0;
        else numOfOutputs = topology[layerCounter+1];
        for (unsigned int neuronCounter = 0; neuronCounter <= topology[layerCounter]; neuronCounter++){
            layers.back().push_back(Neuron(numOfOutputs, neuronCounter));
            std::cout << "made a neuron!"<<"\n";
        }
        layers.back().back().setOutputVal(1.0);
    }

}

void Net::feedForward(const std::vector<double> &inputVals){
    for (unsigned int i = 0; i<inputVals.size(); i++){
        layers[0][i].setOutputVal(inputVals[i]);
    }
    for (unsigned int layerCounter = 1; layerCounter < layers.size(); layerCounter++){
        Layer &prevLayer = layers[layerCounter-1];
        for (unsigned int neuronCounter = 0; neuronCounter<layers[layerCounter].size()-1; neuronCounter++){
            layers[layerCounter][neuronCounter].neuronFeedForward(prevLayer);
        }
    }
}

void Net::backProp(const std::vector<double> &targetVals){
    
    //error calculation
    Layer &outputLayer = layers.back();
    error = 0.0;

    for (unsigned int n = 0; n<outputLayer.size() - 1; ++n){
        double delta = targetVals[n] - outputLayer[n].getOutputVal();
        error+=delta*delta;
    }
    error /= outputLayer.size() -1;
    error = sqrt(error);

    recentAverageError = (recentAverageError * recentAverageSmoothigFactor + error)
                        /(recentAverageSmoothigFactor +1 );

    //output gradient
    for (unsigned int n = 0; n<outputLayer.size() -1; n++){
        outputLayer[n].calcOutputGradient(targetVals[n]);
    }

    //hiddn gradients
    for (unsigned int layerCounter = layers.size()-2; layerCounter>0; --layerCounter){
        Layer &hiddenLayer = layers[layerCounter];
        Layer &nextLayer = layers[layerCounter+1];

        for (unsigned int neuronCounter = 0; neuronCounter<hiddenLayer.size();++neuronCounter){
            hiddenLayer[neuronCounter].calcHiddenGradient(nextLayer);
        }
    }

    //update connection weights
    for (unsigned int layerCounter = layers.size() - 1; layerCounter>0; --layerCounter){
        Layer &currLayer = layers[layerCounter];
        Layer &prevLayer = layers[layerCounter-1];

        for (unsigned int neuronCounter = 0; neuronCounter<currLayer.size() -1 ;neuronCounter++){
            currLayer[neuronCounter].updateInputWeights(prevLayer);
        }
    }
}

//main-------------------------------------------------------------------------

int main(int argc, char *argv[]){

    std::vector<unsigned> topology;
    topology.push_back(4);
    topology.push_back(16);
    topology.push_back(4);
    std::vector<double> inputVals;
    std::vector<double> targetVals;
    std::vector<double> resultVals;
    Net myNet(topology);

    int trainingCounter = 0;
    double val;
    std::ifstream dataFile;
    std::ofstream outFile;

    dataFile.open("training_data.txt");
    outFile.open("NeuralOut.txt");

    while (!dataFile.eof()){
        inputVals.clear();
        targetVals.clear();
        trainingCounter++;
        outFile << "Training #" << trainingCounter << "\n";
        for (unsigned int i=0; i<4; i++){
            dataFile >> val;
            inputVals.push_back(val);
        }
        for (unsigned int i=0; i<4; i++){
        dataFile >> val;
        targetVals.push_back(val);
        }
        myNet.feedForward(inputVals);
        outFile << "Inputs:" << inputVals.at(0) << " "
                << inputVals.at(1)<< " " << inputVals.at(2)<< 
                " " << inputVals.at(3)<< "\n";
        myNet.getResults(resultVals);
        outFile << "Outputs:" << resultVals.at(0) << " "
                << resultVals.at(1) << " " << resultVals.at(2) << " "
                << resultVals.at(3) <<"\n";
        myNet.backProp(targetVals);
        outFile << "Expected Outputs:" << targetVals.at(0) << " " 
                << targetVals.at(1) << " " << targetVals.at(2) << " "  
                << targetVals.at(3) << "\n";

    }

}