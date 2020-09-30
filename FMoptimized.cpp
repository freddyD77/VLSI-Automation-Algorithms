
// C++ program to demonstrate  
// accessing of data members  

#include <iostream> 
#include <fstream>
#include <sstream>
#include <string>
#include<vector>
#include <unordered_map>
#include <unordered_set>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>
#include <chrono> 
using namespace std::chrono;
using namespace std;
//using std::vector;

class Net;
class netList;

class Cell
{
    // Access specifier 
public:

    // Data Members 
    int gain;
    int block;
    bool locked;
    int size;
    //netList *netlist;
    vector<Net*> netlist;


    Cell(int block1, int size1);
    void addNet(Net* inNet);
    void updateGain();

    void initGain() {
        gain = -1 * netlist.size();
    }

    void lock() {
        locked = true;
    }

    void unlock() {
        locked = false;
    }

    void switchBlock() {
        if (block == 0)
            block = 1;
        else if (block == 1)
            block = 0;
    }

    void incGain() {
        gain++;
    }

    void decGain() {
        gain--;
    }

};

class gainTable;

class Net {
public:
    vector<Cell*> B0cellList;
    vector<Cell*> B1cellList;
    //cellList *B0cellList;
    //cellList *B1cellList;

    Net() {
        //B0cellList = new cellList();
        //B1cellList = new cellList();

    }

    void addNode(Cell* cellin) {
        if (cellin->block == 0)
            B0cellList.push_back(cellin);
        else if (cellin->block == 1)
            B1cellList.push_back(cellin);

    }

    void printB0list() {
        for (int x = 0; x < B0cellList.size(); x++) {
            cout << B0cellList.at(x) << " ";
        }
    }

    void updateGains(int previousBlock, Cell* movedCell, gainTable* B0table, gainTable* B1table, int& cutsetSize);


};

Cell::Cell(int block1, int size1) {
    gain = 0;
    block = block1;
    locked = false;
    size = size1;

}

void Cell::addNet(Net* inNet) {
    netlist.push_back(inNet);

}

void Cell::updateGain() {
    int numOfNetsWithEmptyT = 0;
    int numOfNetsWith1F = 0;
    vector<Cell*>* fromList;
    vector<Cell*>* toList;
    for (auto* x : netlist) {
        if (block == 0) {
            fromList = &x->B0cellList;
            toList = &x->B1cellList;
        }
        else {
            toList = &x->B0cellList;
            fromList = &x->B1cellList;
        }
        if (toList->size() == 0)
            numOfNetsWithEmptyT++;
        if (fromList->size() == 1)
            numOfNetsWith1F++;
    }
    gain = numOfNetsWith1F - numOfNetsWithEmptyT;
}

template<typename S>
auto select_random(const S& s, size_t n) {
    auto it = std::begin(s);
    // 'advance' the iterator n times
    std::advance(it, n);
    return it;
}


class gainTable
{
    // Access specifier 
public:

    // Data Members 
    int size;
    int maxGain;
    int currentHighestGain;
    unordered_map<int, unordered_set<Cell*>> table;


    gainTable(int maxgain) {
        maxGain = maxgain;
        currentHighestGain = maxgain;
        size = 0;
        unordered_set<Cell*> empty{};
        for (int x = (-1 * maxgain); x <= maxgain; x++) {
            table[x] = empty;
        }
    }

    void insert(Cell* cell) {
        table[cell->gain].insert(cell);
        size++;
    }

    void insert1(Cell* cell) {
        table[cell->gain].insert(cell);
    }

    Cell* getLargestGainCell() {
        Cell* returnedCell;
        while (1) {
            if (table[currentHighestGain].size() != 0)
                break;
            else if (currentHighestGain == (-1 * maxGain))
                cout << "cannot find a cell in this gainTable\r\n";
            currentHighestGain--;
        }
        auto r = rand() % table[currentHighestGain].size(); // not _really_ random
        returnedCell = *select_random(table[currentHighestGain], r);
        //returnedCell = table[currentHighestGain][0];
        return returnedCell;
    }

    void remove(Cell* cell) {
        //vector<int>::iterator it;
        //it = table[cell->gain].begin();
        //auto i = table[cell->gain].begin();
        table[cell->gain].erase(cell);//assume cell is at position 0 since we always grab the cell from position 0
        size--;
    }

    void updateCellPosition(Cell* cell, int oldgain) {
        table[oldgain].erase(cell);/////hopefully erases cell from this gain bucket (try a print before and after)
        table[cell->gain].insert(cell);
    }


};

void Net::updateGains(int previousBlock, Cell* movedCell, gainTable* B0table, gainTable* B1table, int& cutsetSize) {
    vector<Cell*>* fromList;
    vector<Cell*>* toList;
    gainTable* Btable;
    if (previousBlock == 0) {
        fromList = &B0cellList;
        toList = &B1cellList;
    }
    else {
        fromList = &B1cellList;
        toList = &B0cellList;
    }

    //BEFORE MOVE
    if (toList->size() == 0) {
        for (auto x : *fromList) {
            if (!x->locked) {
                if (x->block == 0)
                    Btable = B0table;
                else
                    Btable = B1table;
                x->incGain();
                Btable->updateCellPosition(x, x->gain - 1);
                if (Btable->currentHighestGain < x->gain)
                    Btable->currentHighestGain = x->gain;
            }
        }
        cutsetSize++;//not sure if using this pointer correctly
    }
    else if (toList->size() == 1) {
        if (!toList->at(0)->locked) {
            if (toList->at(0)->block == 0)
                Btable = B0table;
            else
                Btable = B1table;
            toList->at(0)->decGain();
            Btable->updateCellPosition(toList->at(0), toList->at(0)->gain + 1);
        }
    }

    //ACTUAL MOVE
    fromList->erase(std::remove(fromList->begin(), fromList->end(), movedCell), fromList->end());/////hopefully erases cell from this gain bucket (try a print before and after)
    toList->push_back(movedCell);

    //AFTER MOVE
    if (fromList->size() == 0) {
        for (auto x : *toList) {
            if (!x->locked) {
                if (x->block == 0)
                    Btable = B0table;
                else
                    Btable = B1table;
                x->decGain();
                Btable->updateCellPosition(x, x->gain + 1);
            }
        }
        cutsetSize--;//not sure if using this pointer correctly
    }
    else if (fromList->size() == 1) {
        if (!fromList->at(0)->locked) {
            if (fromList->at(0)->block == 0)
                Btable = B0table;
            else
                Btable = B1table;
            fromList->at(0)->incGain();
            Btable->updateCellPosition(fromList->at(0), fromList->at(0)->gain - 1);
            if (Btable->currentHighestGain < fromList->at(0)->gain)
                Btable->currentHighestGain = fromList->at(0)->gain;
        }
    }
}

void firstPass(gainTable* B0table, gainTable* B1table, unordered_set<Cell*>* lockedCells, float ratio, float totalCells, int& cutsetSize) {
    int previousBlock;
    while (B0table->size > (totalCells * ratio)) {
        Cell* baseCell = B0table->getLargestGainCell();
        previousBlock = baseCell->block;
        baseCell->lock();
        lockedCells->insert(baseCell);
        baseCell->switchBlock();
        B0table->remove(baseCell);
        B1table->size++;
        Net* chosenNet;
        for (int x = 0; x < baseCell->netlist.size(); x++) {
            chosenNet = baseCell->netlist[x];
            chosenNet->updateGains(previousBlock, baseCell, B0table, B1table, cutsetSize);
        }
        
    }
}

void freeLoackedCells(gainTable* B0table, gainTable* B1table, unordered_set<Cell*>* lockedCells) {
    int previousGain;
    gainTable* Btable;
    for (auto* x : *lockedCells) {
        previousGain = x->gain;
        x->unlock();
        x->updateGain();
        if (x->block == 0)
            Btable = B0table;
        else
            Btable = B1table;
        Btable->insert1(x);
    }
    lockedCells->clear();
}

void Pass(gainTable* B0table, gainTable* B1table, unordered_set<Cell*>* lockedCells, float ratio, int& cutsetSize, int& bestCutsetSize, float totalCells) {
    int count = 0;
    int previousBlock;
    gainTable* previousTable;
    gainTable* nextTable;
    Cell* baseCell;

    bestCutsetSize *= 2;
    while (lockedCells->size() < (totalCells * ratio) - 1) {
        count++;
        if (count % 2 == 0) {
            previousTable = B0table;
            nextTable = B1table;
        }
        else {
            previousTable = B1table;
            nextTable = B0table;
        }
        baseCell = previousTable->getLargestGainCell();
        previousBlock = baseCell->block;
        baseCell->lock();
        lockedCells->insert(baseCell);
        baseCell->switchBlock();
        previousTable->remove(baseCell);
        nextTable->size++;
        for (auto* x : baseCell->netlist)
            x->updateGains(previousBlock, baseCell, B0table, B1table, cutsetSize);
        if (cutsetSize < bestCutsetSize)
            bestCutsetSize = cutsetSize;
    }
}


void FMalgoOptimized(float ratio, string sizefile, string netfile, ofstream* outfile, int iter1) {

    auto start = high_resolution_clock::now();
    srand(time(NULL));

    
    /**/
    unordered_map<int, Cell*> cellmap;
    unordered_map<int, Cell*> padmap;
    //float ratio = 0.5;
    //string sizefile = "C:\\temp\\ibm18.are";
    //string netfile = "C:\\temp\\ibm18.net";

    //CREATE ALL CELLS
    //cout << "create new cells\r\n";
    string name;
    string firstLetter;
    int cellsize;
    int cellnum = 0;
    ifstream infile;
    infile.open(sizefile);
    if (!infile)
        cout << "could not open file\r\n";
    while (infile >> name >> cellsize) {
        //stringstream the_name(name.erase(0));
        firstLetter = name.substr(0, 1);
        name.erase(0, 1);
        cellnum = stoi(name);
        //cout << name << " " << firstLetter << " " << cellnum << "\r\n";
        if (firstLetter._Equal("a"))
            cellmap.insert(make_pair(cellnum, new Cell(0, cellsize)));
        else
            padmap.insert(make_pair(cellnum, new Cell(0, cellsize)));
    }
    infile.close();

    //CREATE ALL NETS AND PLUG NET TO CELLS AND CELLS TO NETS
    //cout << "create new nets\r\n";
    int count = 0;
    string startNet;
    Cell* chosenCell;
    Net* currentNet = nullptr;
    unordered_map<int, Cell*>* the_map;

    infile.open(netfile);
    if (!infile)
        cout << "could not open file\r\n";
    std::string line;
    while (std::getline(infile, line))
    {
        if (count >= 5) {
            std::istringstream iss(line);
            if (!(iss >> name >> startNet)) { break; } // error
            firstLetter = name.substr(0, 1);
            name.erase(0, 1);
            cellnum = stoi(name);
            if (startNet._Equal("s"))
                currentNet = new Net();
            if (firstLetter._Equal("a"))
                the_map = &cellmap;
            else
                the_map = &padmap;
            chosenCell = the_map->at(cellnum);
            chosenCell->addNet(currentNet);
            currentNet->addNode(chosenCell);
            //cout << name << " " << startNet << "\r\n";
        }
        else
            count++;
    }
    infile.close();

    /*count = 0;
    unordered_set<Net*> the_nets;
    for (auto x : cellmap) {
        //cout << x.first << " " << x.second << endl;
        for (auto y : x.second->netlist)
            the_nets.insert(y);
    }
    cout << "count is " << the_nets.size() << "\r\n";*/

    //UPDATE GAINS FOR ALL CELLS BASED OFF OF BEING IN THE SAME PARTITION, 
    //AND FIND LARGEST AMOUNT OF NETS FOR A CELL, HENCE THE HIGHEST AVAILABLE GAIN
    //cout << "find largest gain\r\n";
    int largestGain = 0;
    int size = 0;
    for (auto x : cellmap) {
        chosenCell = x.second;
        chosenCell->initGain();
        size = -1 * chosenCell->gain;
        if (size > largestGain)
            largestGain = size;
    }
    for (auto x : padmap) {
        chosenCell = x.second;
        chosenCell->initGain();
    }
    //cout << "largest gain is " << largestGain << "\r\n";

    //CREATE HASH TABLES AND INSERT ALL NODES INTO B0
    //cout << "create the hashtables\r\n";
    gainTable* B0table = new gainTable(largestGain);
    gainTable* B1table = new gainTable(largestGain);
    B0table->currentHighestGain = largestGain;
    B1table->currentHighestGain = largestGain;
    for (auto x : cellmap) {
        B0table->insert(x.second);
    }
    for (auto x : padmap) {
        B0table->insert(x.second);
    }
    float totalCells = B0table->size;
    //cout << totalCells;

    //PERFORM FIRST PASS TO MAKE THE FIRST PARTITION
    int cutsetSize = 0;
    unordered_set<Cell*> lockedCells{};
    //cout << "begin first pass\r\n";
    firstPass(B0table, B1table, &lockedCells, ratio, totalCells, cutsetSize);
    //cout << "finished first pass!\r\n";
    //cout << "b0 size is" << B0table->size << "\r\n";
    //cout << "b1 size is" << B1table->size << "\r\n";
    //cout << "cutset size is" << cutsetSize << "\r\n";
    //cout << "lockedCells size is" << lockedCells.size() << "\r\n";
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    //cout << "end time is " << duration.count() << " milliseconds" << "\r\n";

    //FREE LOCKED CELLS AND UPDATE POSITIONING FOR ALL LOCKED CELLS
    freeLoackedCells(B0table, B1table, &lockedCells);
    //cout << "lockedCells size is" << lockedCells.size() << "\r\n";

    //PERFORM PASSES UNTILL CUTSETSIZE NO LONGER IMPROVES
    count = 0;
    int bestCutesetSize = cutsetSize;
    int previousbestCutsetSize = bestCutesetSize;
    while (bestCutesetSize <= previousbestCutsetSize) {
        previousbestCutsetSize = bestCutesetSize;
        Pass(B0table, B1table, &lockedCells, ratio, cutsetSize, bestCutesetSize, totalCells);
        //cout << "best cutsetsize for this pass is " << bestCutesetSize << "\r\n";
        freeLoackedCells(B0table, B1table, &lockedCells);
        count++;
    }
    //cout << "best cutsetsize is " << previousbestCutsetSize << "\r\n";
    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>(stop2 - start);
    //cout << "end time is " << duration2.count() << " milliseconds" << "\r\n";
    cout << "ratio is " << ratio << " file is " << netfile << " the best cutsetSize is " << bestCutesetSize << " number of passes is "
        << count << " execution time(ms) is " << duration2.count() << "\r\n";
    *outfile << to_string(ratio) + ", " + to_string(iter1) + ", " + to_string(bestCutesetSize) + ", " + to_string(count) + ", " + to_string(duration2.count()) + "\n";

}

int main() {

    

    std::ofstream outfile;
    outfile.open("C:\\temp\\FMalgoOPTdata3.csv");

    float ratio[] = { 0.5, 0.45, 0.4, 0.35 };
    string sizefile[] = { "C:\\temp\\ibm01.are", "C:\\temp\\ibm02.are", "C:\\temp\\ibm03.are", "C:\\temp\\ibm04.are", "C:\\temp\\ibm05.are", "C:\\temp\\ibm06.are"
                        , "C:\\temp\\ibm07.are" , "C:\\temp\\ibm08.are" , "C:\\temp\\ibm09.are" , "C:\\temp\\ibm10.are" , "C:\\temp\\ibm11.are" , "C:\\temp\\ibm12.are"
                        , "C:\\temp\\ibm13.are" , "C:\\temp\\ibm14.are" , "C:\\temp\\ibm15.are" , "C:\\temp\\ibm16.are" , "C:\\temp\\ibm17.are" , "C:\\temp\\ibm18.are" };
    string netfile[] = { "C:\\temp\\ibm01.net", "C:\\temp\\ibm02.net", "C:\\temp\\ibm03.net", "C:\\temp\\ibm04.net", "C:\\temp\\ibm05.net", "C:\\temp\\ibm06.net"
                        , "C:\\temp\\ibm07.net" , "C:\\temp\\ibm08.net" , "C:\\temp\\ibm09.net" , "C:\\temp\\ibm10.net" , "C:\\temp\\ibm11.net" , "C:\\temp\\ibm12.net"
                        , "C:\\temp\\ibm13.net" , "C:\\temp\\ibm14.net" , "C:\\temp\\ibm15.net" , "C:\\temp\\ibm16.net" , "C:\\temp\\ibm17.net" , "C:\\temp\\ibm18.net" };

    for (int y = 16; y < 18/*sizeof(sizefile)/sizeof(sizefile[0])*/; y++) {
        for (int x = 0; x < sizeof(ratio) / sizeof(ratio[0]); x++) {
            FMalgoOptimized(ratio[x], sizefile[y], netfile[y], &outfile, y);
            //cout << x << " " << y << "\r\n";
        }
    }

    outfile.close();



    return 0;

    return 0;
}