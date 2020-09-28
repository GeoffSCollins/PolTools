#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <unordered_map>

using namespace std;

vector<string> splitString(const string& str, const string& delim) {
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

string stepOneAndTwo(string bedFilename, int maxFragmentSize) {
    // Output is chromosome, left, right, lenSum, sum, strand

    map<string, map<string, map<string, unordered_map<int, int>>>> countdata;

    ifstream bedFile(bedFilename);

    string chromosome, strLeft, strRight, name, score, strand;
    int left, right;

    // The strand is not getting updated for some reason!!!!
    while (bedFile >> chromosome >> left >> right >> name >> score >> strand) {
        int fivePrimeEnd = 0;

        if (strand.compare("+") == 0) {
            fivePrimeEnd = left;
        } else {
            fivePrimeEnd = right-1;
        }

        // Check if it is in the dictionary
        if (countdata[chromosome]["five"][strand].find(fivePrimeEnd) == countdata[chromosome]["five"][strand].end()) {
            // We need to add the key
            countdata[chromosome]["five"][strand][fivePrimeEnd] = 0;
        }

        if (countdata[chromosome]["length"][strand].find(fivePrimeEnd) == countdata[chromosome]["length"][strand].end()) {
            // We need to add the key
            countdata[chromosome]["length"][strand][fivePrimeEnd] = 0;
        }

        countdata[chromosome]["five"][strand][fivePrimeEnd]++;

        if ((right - left) <= maxFragmentSize) {
            countdata[chromosome]["length"][strand][fivePrimeEnd] += (right - left);
        }
    }

    bedFile.close();

    // Output the data
    string outputFilename = bedFilename.replace(bedFilename.find(".bed"), sizeof(".bed") - 1, "-2-output.txt");
    ofstream outputFile(outputFilename);

    // The countdata map is now filled, so we output all of the data
    for(auto const& chromMap: countdata) {
        string chrom = chromMap.first;

        for(auto const& strandMap: countdata[chrom]["five"]) {
            string strand = strandMap.first;

            for(auto const& positionMap: countdata[chrom]["five"][strand]) {
                int position = positionMap.first;

                // Only print out positions which have data
                if (countdata[chrom]["five"][strand].find(position) != countdata[chrom]["five"][strand].end()) {
                    outputFile << chrom << "\t" << position << "\t" << position + 1 << "\t" << countdata[chrom]["length"][strand][position];
                    outputFile << "\t" << countdata[chrom]["five"][strand][position] << "\t" << strand << endl;
                }
            }
        }
    }

    outputFile.close();

    return outputFilename;
}


map<string, vector<string>> readStepTwoFile(string fileName, string &strand) {
    ifstream file(fileName);
    string chrom, left, right, readLengthSum, readCount;
    string line;
    map<string, vector<string>> dataByChromosome;

    while (file >> chrom >> left >> right >> readLengthSum >> readCount >> strand) {
        if (dataByChromosome.find(chrom) == dataByChromosome.end()) {
            // The key does not exist
            vector<string> regions;
            dataByChromosome.insert({chrom, regions});
        }

        // Add the current region to the map
        dataByChromosome[chrom].push_back(left + "\t" + right + "\t" + readLengthSum + "\t" + readCount + "\t" + strand);
    }

    file.close();
    return dataByChromosome;
}


void buildMap(unordered_map<string, vector<int>> &tsrMap, vector<string> regions, int stepSize) {
    // For each region in the map values (given from step 2 file)
    for (string region: regions) {
        vector<string> splitRegion = splitString(region, "\t");
        // Split the string by tabs to grab the entries

        int tss = stoi(splitRegion[0]);
        int tssPlusOne = stoi(splitRegion[1]);
        int readLenSum = stoi(splitRegion[2]);
        int readCount = stoi(splitRegion[3]);

         // For each sub region of length stepsize starting at the tss + 1 - the stepSize to the tss + 1
         // This is the window loop?
         for (int position = tss + 1 - stepSize; position < tss + 1; position++) {
            int startRegion = position;
            int endRegion = position + stepSize;

            // Don't deal with negative positions
            if (startRegion < 0) {
                continue;
            }

            string mapKey = to_string(startRegion) + "-" + to_string(endRegion);

            if (tsrMap.find(mapKey) == tsrMap.end()) {
                // This means the TSR was not added to the map
                vector<int> mapData = {0, 0, 0, 0, 0, 0};
                tsrMap[mapKey] = mapData;
            }

            tsrMap[mapKey][0] += readLenSum;
            tsrMap[mapKey][1] += readCount;

            if (readCount >= tsrMap[mapKey][3]) {
                tsrMap[mapKey][2] = tss;
                tsrMap[mapKey][3] = readCount;
            }

            // This is a helper for the avgTSS
            tsrMap[mapKey][4] += tssPlusOne * readCount;

            // This is a helper for the stdev
            tsrMap[mapKey][5] += tssPlusOne * readCount;
         }
    }
}


string stepThree(string stepTwoFilename, int stepSize, int minSeqDepth, int minAvgTranscriptLength) {
    string strand;

    map<string, vector<string>> dataByChromosome = readStepTwoFile(stepTwoFilename, strand);

    map<string, map<int, string>> sortedOutputData;

    // Loop over keys
    for(auto const& key: dataByChromosome) {
        string chrom = key.first;
        vector<string> regions = key.second;
        unordered_map<string, vector<int>> tsrMap;

        buildMap(tsrMap, regions, stepSize);

        for (auto const& subMap: tsrMap) {
            string subMapKey = subMap.first;

            int windowStart = stoi(splitString(subMapKey, "-")[0]);
            int windowEnd = stoi(splitString(subMapKey, "-")[1]);

            int readLenSum = tsrMap[subMapKey][0];
            int readCount = tsrMap[subMapKey][1];

            int maxTSS = tsrMap[subMapKey][2];
            int maxReadCount = tsrMap[subMapKey][3];
            double avgTSSHelper = tsrMap[subMapKey][4];

            // Only keep the ones with RCOVSUM >= MINSEQDEPTH and (RLSUM / RCOVSUM) >= AVGTRANSLEN
            if (readCount >= minSeqDepth && (readLenSum / readCount) >= minAvgTranscriptLength) {
                // We need to get the other values and then print it out
                double averageTSS = avgTSSHelper / readCount;

                string data = chrom + "\t" + to_string(windowStart) + "\t" + to_string(windowEnd) + "\t" + to_string(readLenSum) + "\t" + to_string(readCount);
                data += "\t" + strand + "\t" + to_string(maxTSS) + "\t" + to_string(maxTSS + 1) + "\t" + to_string(maxReadCount) + "\t";
                data += to_string(averageTSS) + "\t" + to_string(maxTSS + 1 - averageTSS);

                sortedOutputData[chrom][windowStart] = data;
            }
        }
    }

    // Sort the helper file by chromosome and then by left position
    string outputFilename = stepTwoFilename.replace(stepTwoFilename.find("-2-output.txt"), sizeof("-2-output.txt") - 1, "-3-output.txt");
    ofstream outputFile(outputFilename);

    for (auto const& chromMap: sortedOutputData) {
        string chrom = chromMap.first;

        for (auto const& positionMap: sortedOutputData[chrom]) {
            int position = positionMap.first;
            outputFile << sortedOutputData[chrom][position] << endl;
        }
    }


    outputFile.close();
    return outputFilename;
}

int main(int argc, char** argv)
{
    if (argc != 7) {
        return 1;
    }

    string bedFilename = argv[1];
    int windowSize = stoi(argv[2]);
    int minSeqDepth = stoi(argv[3]);
    int minAvgTranscriptLength = stoi(argv[4]);
    int maxFragmentSize = stoi(argv[5]);
    string chromSizeFilename = argv[6];


    string stepTwoOutputFilename = stepOneAndTwo(bedFilename, maxFragmentSize);
    stepThree(stepTwoOutputFilename, windowSize, minSeqDepth, minAvgTranscriptLength);

    return 0;
}

