#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <array>
#include <vector>
#include <list>
#include <algorithm>
#include <queue>

#define DB_MAX_SIZE 100


using namespace std;

// Read from library file and record each gate's cell_delay and output_slew table in myMap
map<string, array<array<array<double, 7>, 7>, 2> > myMap;
map<string, double> myCap; // store input capacitancec for each gate
map<string, array<array<double, 7>, 2> > myIndex;
map<int, string> intToStringMap; // Recorder to find each gate in myMap
map<int, double> arriveTimes; // record arrival time for each gate
map<int, double> slewValue; // record slew value for each gate
map<int, double> maxDelayValues; // record min delay value for each gate *****

double maxCircuitDelay = 0;
vector<int> criticalPath;
// Map to store the required arrival time for each gate
map<int, double> requiredArrivalTimes;
// Map to store the slack for each gate
map<int, double> slackTime;
int minSlackPrimaryOutputGateID = -1;


int counter = 1;

// Read from circuit file and initialize size to 1000
vector< list<int> > gates(1000);
vector< list<int> > fanouts(1000);


int parseLib(char *fName);
int parseCircuit(char *fName);

struct ParenCommaEq_is_space : std::ctype<char> {
    ParenCommaEq_is_space() : std::ctype<char>(get_table()) {}
    static mask const* get_table()
    {
        static mask rc[table_size];
        rc['('] = std::ctype_base::space;
        rc[')'] = std::ctype_base::space;
        rc[','] = std::ctype_base::space;
        rc['='] = std::ctype_base::space;
        rc[' '] = std::ctype_base::space;
        rc['\t'] = std::ctype_base::space;
        return &rc[0];
        }
};
// Do delay value interpolation
double delayInterpolate(string gateType, double inputSlew, double outLoad){

    int i, j;
    array<double, 7> tmpSlewTable = myIndex[gateType][0];
    array<double, 7> tmpLoadTable = myIndex[gateType][1];
    for (i = 0; i < 6; i++) {
        if (tmpSlewTable[i] <= inputSlew && inputSlew <= tmpSlewTable[i+1]) break;
    }
    
    for (j = 0; j < 6; j++) {
        if (tmpLoadTable[j] <= outLoad && outLoad <= tmpLoadTable[j+1]) break;
    }

    double t1 = tmpSlewTable[i];
    double t2 = tmpSlewTable[i+1];
    double C1 = tmpLoadTable[j];
    double C2 = tmpLoadTable[j+1];
    
    array<array<double, 7>, 7> tmpDelayTable = myMap[gateType][0];
    double v11 = tmpDelayTable[i][j];
    double v12 = tmpDelayTable[i][j+1];
    double v21 = tmpDelayTable[i+1][j];
    double v22 = tmpDelayTable[i+1][j+1];

    double value = v11*(C2 - outLoad)*(t2 - inputSlew) + v12*(outLoad - C1)*(t2 - inputSlew) 
                  + v21*(C2 - outLoad)*(inputSlew - t1) + v22*(outLoad - C1)*(inputSlew - t1);
    value /= (C2 - C1)*(t2 - t1);
    return value;

}
// Do slew value interpolation
double slewInterpolate(string gateType, double inputSlew, double outLoad){

    int i, j;
    array<double, 7> tmpSlewTable = myIndex[gateType][0];
    array<double, 7> tmpLoadTable = myIndex[gateType][1];
    for (i = 0; i < 6; i++) {
        if (tmpSlewTable[i] <= inputSlew && inputSlew <= tmpSlewTable[i+1]) break;
    }
    
    for (j = 0; j < 6; j++) {
        if (tmpLoadTable[j] <= outLoad && outLoad <= tmpLoadTable[j+1]) break;
    }

    double t1 = tmpSlewTable[i];
    double t2 = tmpSlewTable[i+1];
    double C1 = tmpLoadTable[j];
    double C2 = tmpLoadTable[j+1];
    
    array<array<double, 7>, 7> tmpSlewValueTable = myMap[gateType][1];
    double v11 = tmpSlewValueTable[i][j];
    double v12 = tmpSlewValueTable[i][j+1];
    double v21 = tmpSlewValueTable[i+1][j];
    double v22 = tmpSlewValueTable[i+1][j+1];

    double value = v11*(C2 - outLoad)*(t2 - inputSlew) + v12*(outLoad - C1)*(t2 - inputSlew) 
                  + v21*(C2 - outLoad)*(inputSlew - t1) + v22*(outLoad - C1)*(inputSlew - t1);
    value /= (C2 - C1)*(t2 - t1);

    return value;

}
// Compute each gate delay and slew value
void computeGateDelay(int gateIndex) {

    string gateType = intToStringMap[gates[gateIndex].front()];

    double maxArrive = 0.0;
    double maxSlewValue = 0.0;
    int count = gates[gateIndex].size()-1; // record # of inputs

    // Calculate output load, the same for each input
    double outLoad = 0.0;

    // If fanouts[input] is empty, then it is an output
    if (fanouts[gateIndex].empty()) {
        outLoad = 4 * myCap["INV"];
    } else {
        // Otherwise, sum up the capacitances of all fanout gates
        for (int output : fanouts[gateIndex]) {
            string outType = intToStringMap[gates[output].front()];
            if (outType == "DFF"){
                outLoad += 4 * myCap["INV"];
            }else{
                outLoad += myCap[outType];
            }
        }
    }
    int y = -1;
    double firstLoad = myIndex[gateType][1][0];
    double endLoad = myIndex[gateType][1][6];
    if (outLoad < firstLoad){
        y = 0;
    }else if(outLoad > endLoad){
        y = 6;
    }else{
        // Check for outLoad in myIndex[1]
        for (int i = 0; i < 7; ++i) {
            if (myIndex[gateType][1][i] == outLoad) {
                y = i;
                break;
            }
        }
    }

    // Calculate delay for each input

    auto it = gates[gateIndex].begin();
    ++it; // Skip the gate type
    for (; it != gates[gateIndex].end(); ++it) {
        int input = *it;

        // input slew
        double inputSlew = 0.0;
        
        // Check if it's a primary input (size of the linked list is 1)
        if (gates[input].size() == 1) {
            inputSlew = 0.002; // ps to ns
        } else {
            // Otherwise, fetch the slewValue for that gate
            inputSlew = slewValue[input];
        }


        // calculate x and y based on inputSlew and outLoad
        int x = -1;
        double curDelay = 0;
        double curSlewValue = 0;

        // Check the edge cases
        double firstSlew = myIndex[gateType][0][0];
        double endSlew = myIndex[gateType][0][6];

        if (inputSlew < firstSlew){
            x = 0;
        }else if(inputSlew > endSlew){
            x = 6;
        }else{
            // Check for inputSlew in myIndex[0]
            for (int i = 0; i < 7; i++) {
                if (myIndex[gateType][0][i] == inputSlew) {
                    x = i;
                    break;
                }
            }   
        }
        // Four different interpolation cases
        if (x == -1 and y == -1) {
            curDelay = delayInterpolate(gateType, inputSlew, outLoad);
            curSlewValue = slewInterpolate(gateType, inputSlew, outLoad);
        }else if(x == -1){
            curDelay = delayInterpolate(gateType, inputSlew, myIndex[gateType][1][y]);
            curSlewValue = slewInterpolate(gateType, inputSlew, myIndex[gateType][1][y]);
        }else if(y == -1){
            curDelay = delayInterpolate(gateType, myIndex[gateType][0][x], outLoad);
            curSlewValue = slewInterpolate(gateType, myIndex[gateType][0][x], outLoad);
        }else{
            curDelay = myMap[gateType][0][x][y];
            curSlewValue = myMap[gateType][1][x][y];
        }

        // calculate Arrive time and update longest one
        double arriveTime = 0.0;
        if (gates[input].size() != 1){
            arriveTime = arriveTimes[input];
        }else{
            arriveTimes[input] = 0.0;
        }
        if (count >= 2){
            curDelay = curDelay * (count/2.0);
            curSlewValue =  curSlewValue * (count/2.0);
            arriveTime += curDelay;
        }else{
            arriveTime += curDelay;
        }

        if (arriveTime > maxArrive){
            maxArrive = arriveTime;
            maxSlewValue = curSlewValue;
            maxDelayValues[gateIndex] = curDelay;
        }
    }
    //store the gate's arrive time and slew value
    arriveTimes[gateIndex] = maxArrive;
    slewValue[gateIndex] = maxSlewValue ;

    // if it's ouput, update longest whole circuit delay
    if (fanouts[gateIndex].empty()){
        if (maxArrive > maxCircuitDelay){
            minSlackPrimaryOutputGateID = gateIndex;
        }
        maxCircuitDelay = max(maxCircuitDelay, maxArrive);
        requiredArrivalTimes[gateIndex] = maxCircuitDelay; // The way to know which gates are outputs
    }
    
}
// Break DFF into an output pin and an input pin with no connections between them
void breakDFF(){
    // Iterate over all the gates
    for (size_t i = 0; i < gates.size(); ++i) {
        // Check if the gate is a "DFF"
        if (!gates[i].empty() && intToStringMap[gates[i].front()] == "DFF") {
            // Set the first element to -1
            gates[i].front() = -1;
            // Remove all other elements from the list
            gates[i].erase(++gates[i].begin(), gates[i].end());
        }
    }
}
// Forward Traverse (Max circuit Delay)
void topologicalTraversal() {
    int gateSize = gates.size();
    vector<int> inDegree(gateSize, 0);  // let all degree to 0
    queue<int> Q;

    // calculate each gate degree
    for (int i = 0; i < gateSize; i++) {
        if (!gates[i].empty()) { // first value means gate type
            inDegree[i] = gates[i].size()-1; 
            if (inDegree[i] == 0){  // put those degree = 0 in Queue
                Q.push(i);
            }
        }
    }

    while (!Q.empty()) {
        int currentGate = Q.front();
        Q.pop();

        // compute current Gate delay
        string gateType = intToStringMap[gates[currentGate].front()];
        if (gates[currentGate].size() != 1){ // if gateIndex is input, it's not a gate, so just ignore it
                computeGateDelay(currentGate);
        }

        // update degree
        for (int outputGate : fanouts[currentGate]) {
            inDegree[outputGate]--;
            if (inDegree[outputGate] == 0) {
                Q.push(outputGate);
            }
        }
    }

    for (int i = 0; i < gateSize; ++i) {
        if (inDegree[i] != 0) {
            cout << "There's a cycle in the circuit!" << endl;
            return;
        }
    }
}
// Back Traverse (slack time for each gate)
void backTraversal(){
    int sz = gates.size();
    // Step 1: Compute in-degrees for all gates (considering fanouts as fanins for this step)
    vector<int> backInDegrees(sz, 0);
    queue<int> Q;


    // Step 1: Compute in-degrees for all gates (considering fanouts as fanins for this step)
    for (int i = 0; i < sz; i++) {
        backInDegrees[i] = fanouts[i].size();
        if (requiredArrivalTimes[i] != 0) { 
            // To find all primary outputs
            Q.push(i);
        }
    }

    // Step 2: Perform the backward traversal
    while (!Q.empty()) {
        int currentGate = Q.front();
        Q.pop();

        // Process each fanin of the current gate (from the gates list)
        auto it = gates[currentGate].begin();
        ++it;  // Skip the first element which denotes gate type

        for (; it != gates[currentGate].end(); ++it) {
            int inputGate = *it;

            double curRT = requiredArrivalTimes[currentGate] - maxDelayValues[currentGate];
            // Update the required arrival time for the input gate by taking the minimum
            if (requiredArrivalTimes[inputGate] == 0 || curRT < requiredArrivalTimes[inputGate]) {
                requiredArrivalTimes[inputGate] = curRT;
            }

            // Decrement the backInDegrees count
            backInDegrees[inputGate]--;
            if (backInDegrees[inputGate] == 0) {
                Q.push(inputGate);  // Push the inputGate to the queue for further traversal
            }
        }
    }

}
// Compute slack of each gate
void computeSlack(){
    // Traverse through all the gates to compute the slack
    for (const auto& entry : arriveTimes) {
        int id = entry.first;
        if(requiredArrivalTimes[id] != 0.0){
            slackTime[id] = requiredArrivalTimes[id] - entry.second;
        }
    }
};
// Function to find the critical path
vector<int> findCriticalPath(map<int, double>& slackTime, vector< list<int> >& gates) {
    vector<int> path;
    int minSlackGateID = minSlackPrimaryOutputGateID;
    // Add the primary output with minimum slack to the critical path
    path.push_back(minSlackGateID);

    while (gates[minSlackGateID].size() > 1) {
        // Get all inputs for the current gate, starting from the second position
        auto it = gates[minSlackGateID].begin();
        ++it; // Move iterator to the second position
        list<int> inputs(it, gates[minSlackGateID].end());

        // Find the input with the minimum slack
        int minSlackInputID = *min_element(inputs.begin(), inputs.end(),
                                           [&slackTime](int i1, int i2) {
                                               return slackTime[i1] < slackTime[i2];
                                           });

        // Add the input with minimum slack to the critical path
        path.push_back(minSlackInputID);

        // Move to the next gate
        minSlackGateID = minSlackInputID;
    }

    return path;
}
// Output information to kt_traversal.txt
int outputToFile(){
    // Open the file for writing
    ofstream file("ckt_traversal.txt");
    if (!file) {
        cerr << "Unable to open file for writing." << endl;
        return 1;
    }

    // Writing to the file
    file << "Circuit delay: " << maxCircuitDelay*1000 << "ps\n\n";
    file << "Gate slacks:\n";
    
    for (const auto& entry : slackTime) {
        int gateID = entry.first;
        double slack = entry.second;
        if (slack != 0) {
            int gateType = gates[gateID].front();
            string gateTypeName = intToStringMap[gateType];
            file << gateTypeName << "-n" << gateID+1 << ": " << slack*1000 << "ps\n";
        }
    }

    file << "\nCritical path:\n";
    for (auto it = criticalPath.rbegin(); it != criticalPath.rend(); ++it) {
        int gateID = *it;
        int gateType = gates[gateID].front();
        string gateTypeName = intToStringMap[gateType];
        file << gateTypeName << "-n" << gateID+1;
        if (it != criticalPath.rend() - 1) { // if it is the first element
            file << ",";
        }
    }
    file << "\n";

    // Close the file
    file.close();

    cout << "Data written to ckt_traversal.txt successfully." << endl;
    return 0;
}
// Main
int main(int argc, char* argv[]) {
    if (argc == 1) {
        cout << "Please input the LIB file and CIRCUIT file\n";
        return -1;
    }

    if (argc == 2) {
        cout << "Please input the CIRCUIT file\n";
        return -1;
    }

    if (argc > 3) {
        cout << "Only need to input the LIB file and CIRCUIT file\n";
        return -1;
    }
    //First step: Read LIB file and store cell_delay and output_slew Table for each gate in myMap
    
    parseLib(argv[1]);
    parseCircuit(argv[2]);
    intToStringMap[-1] = "INPUT";
    breakDFF(); //Break DFF into an output pin and an input pin with no connections between them
    fanouts.resize(gates.size(), list<int>()); // update size
    for (int i = 0; i < gates.size(); i++) { // update all fanouts[input1] = linkedlist[output1, output2...]
        if (gates[i].size() > 1) {
            int skipCount = 1;
            int count = 0;

            for (int input : gates[i]) {
                if (count >= skipCount) {
                    fanouts[input].push_back(i);
                }
                count++;
            }
        }
    }

    // step 1: start forward traverse the gates by topological way
    topologicalTraversal();

    //Step 2: Compute the required arrival time for each primary output
    double requiredTimeForOutputs = 1.1 * maxCircuitDelay;
    for(int i = 0; i < requiredArrivalTimes.size(); i++) {
        if (requiredArrivalTimes[i]!=0) { // Assuming a gate with no fanouts is a primary output
            requiredArrivalTimes[i] = requiredTimeForOutputs;
        }
    }
    // Step 3: Backward traversal to compute the required arrival time for each gate
    backTraversal();
    // Step 4: Compute the slack for each gate
    computeSlack();
    // Step 5: Determine the critical path
    criticalPath = findCriticalPath(slackTime, gates);
    // Step 6: Print out in "ckt_traversal.txt"
    outputToFile();

    return 0;
}


// Read Circuit file
int parseCircuit(char *fName){

    ifstream ifs(fName);
	if (!ifs.is_open()) {
        cout << "Error opening file " << fName << endl;
		return -1;
	}

    // Keep reading file while no stream error occurs
    // like end of file (eof) etc.

	while (ifs.good()) {
        
        string lineStr;
        getline(ifs, lineStr);
        istringstream ifs(lineStr);
        // containing one line of text from the file
        ifs.imbue(locale(cin.getloc(),new ParenCommaEq_is_space));
        string firstWord;
        ifs >> firstWord;


         // Check if the input info line
        if (firstWord.compare("INPUT") == 0) {
            string second_word;
            ifs >> second_word;
            // the second word will be a input node number
            // We need to initialize each input node's linked list in vector with initialize "-1", which means INPUT
            int gateID = stoi(second_word);

            if (gateID > gates.size()){ // Here means we need to resize the vector
                int thousandsPlace = gateID / 1000;
                if ((thousandsPlace * 1000) == gateID){ //Determine if gateID is already a multiple of 1000 or not
                    gates.resize(gateID, list<int>());
                }else{  // If not, plus 1000 more
                    gates.resize((thousandsPlace + 1) * 1000, list<int>());

                }
            }
            gates[gateID-1].push_back(-1);
            

        }else if(isdigit(firstWord[0])){ // determine whether the first index of is a digital number or not
            // If it is, here we find the line about gate ID
            int gateID = stoi(firstWord);
            // Same as above
            if (gateID > gates.size()){ // Here means we need to resize the vector
                int thousandsPlace = gateID / 1000;
                if ((thousandsPlace * 1000) == gateID){ //Determine if gateID is already a multiple of 1000 or not
                    gates.resize(gateID, list<int>());
                }else{  // If not, plus 1000 more
                    gates.resize((thousandsPlace + 1) * 1000, list<int>());

                }
            }
            

            // Firstly determine which type of gate
            string second_word;
            ifs >> second_word;
            transform(second_word.begin(), second_word.end(), second_word.begin(), ::toupper);

            // Iterate through the intToStringMap to find if there is a key that matches the value of the string second_word.
            bool isFind = false;
            int keyNumber;
            for (const auto& pair : intToStringMap) {
                if (pair.second == second_word) {
                    // Get key number
                    keyNumber = pair.first;
                    isFind = true;
                    break;
                }
            }
            // If not found, record new element
            if (!isFind) {
                intToStringMap[counter] = second_word;
                gates[gateID-1].push_back(counter); // Give linked list an indicator of gate type
                counter++;
            }else{
                gates[gateID-1].push_back(keyNumber); // Give linked list an indicator of gate type
            }
            
            // Secondly record each input node number in the linked list
            string third_word;
            ifs >> third_word;
            gates[gateID-1].push_back(stoi(third_word)-1);

            //Determine if there is an input gate number after it
            while (ifs >> third_word && isdigit(third_word[0])) {
                // Convert extracted strings to integers
                long number;
                number = stoi(third_word);
                gates[gateID-1].push_back(number-1);
            }

        }
	}
    ifs.close();
    return 0;
};


// Read Library File
int parseLib(char *fName) {
    // ------------------------------------------------------------------
    // Part 1. Parsing Lib file
    // ------------------------------------------------------------------
    // Open file. Check if opened correctly
    // Always check for this to save a lot of pain if an error occurs
	ifstream ifs(fName);
	if (!ifs.is_open()) {
        cout << "Error opening file " << fName << endl;
		return -1;
	}

    size_t db_idx = 0;
    string key;

    // Keep reading file while no stream error occurs
    // like end of file (eof) etc.
	while (ifs.good()) {
        string first_word;
        ifs >> first_word;

        // An info entry must start with "cell" as the first word exactly.
        if (first_word.compare("cell") == 0) {

            string second_word;
            ifs >> second_word;

            // The gate name has the format of "(<gate_name>)".
            // So find the index till the closing bracket.
            // Make sure to ignore the first index which has the opening bracket (hence the 1 and -1)
            size_t delim_pos = second_word.find(")");
            string gate_name = second_word.substr(1, delim_pos-1);

            key = gate_name;

            // Create two 7x7 arrays, 0 for cell_delay, 1 for output_slew
            array<array<array<double, 7>, 7>, 2> newArray;
            for (int k = 0; k < 2; k++) {
                for (int i = 0; i < 7; i++) {
                    for (int j = 0; j < 7; j++) {
                        newArray[k][i][j] = 0.0; // default to set all value to 0.0
                    }
                }
            }
            //Inserting key-value pairs into a map
            myMap[key] = newArray;

            // Deal with capacitance
            string tmp;
            getline(ifs, tmp);
            ifs >> tmp; // Ignore "capacitance"
            ifs >> tmp; // Ignore ":"
            ifs >> tmp; // <value>;
            size_t semiC_pos = tmp.find(";");
            string capStr = tmp.substr(0, semiC_pos);
            double cap = stod(capStr);
            //Inserting key-value pairs into a map
            myCap[key] = cap;

        }
        // The cell delays will start after the word "cell_delay(Timing_7_7)" exactly
        else if (first_word.compare("cell_delay(Timing_7_7)") == 0) {
            // Read 3 lines that contain the rest of above match, index 1 and index 2
            string tmp;
            getline(ifs, tmp); 
            getline(ifs, tmp);  // Go to index_1 line

            // The delays will be between " ". Find the opening ".
            size_t start_idx = tmp.find("\"");
            size_t end_idx = tmp.find("\"", start_idx + 1);
            string data_str = tmp.substr(start_idx + 1, end_idx - start_idx - 1);
        
            istringstream data_stream(data_str);
            for (size_t k = 0; k < 7; k++) {
                double idx1;
                char delim;
                data_stream >> idx1 >> delim;
                myIndex[key][0][k] = idx1; // change ns to ps
            }

             getline(ifs, tmp); // Go to index_2 line

            // The delays will be between " ". Find the opening ".
            start_idx = tmp.find("\"");
            end_idx = tmp.find("\"", start_idx + 1);
            data_str = tmp.substr(start_idx + 1, end_idx - start_idx - 1);
        
            istringstream data_stream2(data_str);
            for (size_t m = 0; m < 7; m++) {
                double idx2;
                char delim;
                data_stream2 >> idx2 >> delim;
                myIndex[key][1][m] = idx2;
            }

            // Now finish store index1 and index2

            // From here on the next 7 lines will contain our delays
            for (size_t i = 0; i < 7; i++) {
                getline(ifs, tmp);

                // The delays will be between " ". Find the opening ".
                size_t start_delim_idx = tmp.find("\"");

                // Find the closing ".
                // The second argument is where we want to start our search
                // Ignore the first match so we don't get the same index again
                size_t end_delim_idx = tmp.find("\"", start_delim_idx + 1);

                // The second arg in substr in no. of characters, not the ending index
                string data_str = tmp.substr(start_delim_idx + 1, end_delim_idx - start_delim_idx - 1);

                // Convert this remaining string to a stream so we can parse our data in doubles
                istringstream data_stream(data_str);
                for (size_t j = 0; j < 7; j++) {
                    double delay;
                    char delim;
                    data_stream >> delay >> delim;
                    myMap[key][0][i][j] = delay;


                }
            }

            // At the end of nested for loop we will have finised reading the 7x7 delays.
            // Increment our database pointer so we can store the next entry
            db_idx++;

            // Sanity check for array size
            if (db_idx >= DB_MAX_SIZE) {
                cout << "ERROR: Insufficient map size for saving data" << endl;
                ifs.close();
                return -1;
            }
        }else if(first_word.compare("output_slew(Timing_7_7)") == 0){
            // Read 3 lines that contain the rest of above match, index 1 and index 2
            string tmp;
            getline(ifs, tmp);
            getline(ifs, tmp);
            getline(ifs, tmp);

            // From here on the next 7 lines will contain our delays
            for (size_t i = 0; i < 7; i++) {
                getline(ifs, tmp);

                // The delays will be between " ". Find the opening ".
                size_t start_delim_idx = tmp.find("\"");

                // Find the closing ".
                // The second argument is where we want to start our search
                // Ignore the first match so we don't get the same index again
                size_t end_delim_idx = tmp.find("\"", start_delim_idx + 1);

                // The second arg in substr in no. of characters, not the ending index
                string data_str = tmp.substr(start_delim_idx + 1, end_delim_idx - start_delim_idx - 1);

                // Convert this remaining string to a stream so we can parse our data in doubles
                istringstream data_stream(data_str);
                for (size_t j = 0; j < 7; j++) {
                    double delay;
                    char delim;
                    data_stream >> delay >> delim;

                    myMap[key][1][i][j] = delay;

                }
            }

            // At the end of nested for loop we will have finised reading the 7x7 delays.
            // Increment our database pointer so we can store the next entry
            db_idx++;

            // Sanity check for array size
            if (db_idx >= DB_MAX_SIZE) {
                cout << "ERROR: Insufficient map size for saving data" << endl;
                ifs.close();
                return -1;
            }
        }
	}
    ifs.close();

	return 0;
}


