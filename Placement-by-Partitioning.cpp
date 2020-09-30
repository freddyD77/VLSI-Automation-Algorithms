
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

int nodecnt = 0;
int hmin = 16;
int B0size = 0;

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
    int xpos;
    int ypos;
    vector<Net*> netlist;
    vector<Net> saved_netlist;


    Cell(int block1, int size1, string name);
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

    Net() {
        //outcell = NULL;
        //traversed = false;
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

Cell::Cell(int block1, int size1, string name) {
    gain = 0;
    block = block1;
    locked = false;
    size = size1;
    xpos = 0;
    ypos = 0;

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
    int numCells;
    int size;
    int maxGain;
    int currentHighestGain;
    unordered_map<int, unordered_set<Cell*>> table;


    gainTable(int maxgain) {
        maxGain = maxgain;
        currentHighestGain = maxgain;
        size = 0;
        numCells = 0;
        unordered_set<Cell*> empty{};
        for (int x = (-1 * maxgain); x <= maxgain; x++) {
            table[x] = empty;
        }
    }


    void insert(Cell* cell) {
        table[cell->gain].insert(cell);
        size+=cell->size;
        numCells++;
    }

    void insert1(Cell* cell) {
        table[cell->gain].insert(cell);
        if (cell->gain > currentHighestGain)
            currentHighestGain = cell->gain;
    }

    Cell* getLargestGainCell() {
        Cell* returnedCell;
        int leftoverCellcnt = 0;
        while (1) {
            if (table[currentHighestGain].size() != 0)
                break;
            else if (currentHighestGain == (-1 * maxGain)) {
                cout << "cannot find a cell in this gainTable\r\n";
                cout << "numCells in table is " << numCells << "\r\n";
                cout << "num of nodes moved is " << nodecnt << "\r\n";
                for (auto x : table) {
                    if (x.second.size() != 0)
                        for(auto y : x.second)
                            leftoverCellcnt++;
                }
                cout << "left over cells are " << leftoverCellcnt << "\r\n";
                cout << "gain pointer is at " << currentHighestGain << "\r\n";
                cout << "max gain is " << maxGain << "\r\n";
            }
            currentHighestGain--;
        }
        auto r = rand() % table[currentHighestGain].size(); // not _really_ random
        returnedCell = *select_random(table[currentHighestGain], r);
        //returnedCell = table[currentHighestGain][0];
        return returnedCell;
    }

    void remove(Cell* cell) {
        table[cell->gain].erase(cell);//assume cell is at position 0 since we always grab the cell from position 0
        size-=cell->size;
        numCells--;
    }

    void updateCellPosition(Cell* cell, int oldgain) {
        table[oldgain].erase(cell);/////hopefully erases cell from this gain bucket (try a print before and after)
        table[cell->gain].insert(cell);
    }


};

int lastMovedCellSize = 0;

void Net::updateGains(int previousBlock, Cell* movedCell, gainTable* B0table, gainTable* B1table, int& cutsetSize) {
    vector<Cell*>* fromList;
    vector<Cell*>* toList;
    gainTable* Btable;
    lastMovedCellSize = movedCell->size;
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
    /*fromList=NULL;
    toList=NULL;
    Btable=NULL;*/
}

void firstPass(gainTable* B0table, gainTable* B1table, unordered_set<Cell*>* lockedCells, float ratio, float totalCells, int& cutsetSize) {
    int previousBlock;
    while (B0table->size > (totalCells * ratio) && B0table->numCells > 1) {//B0 is smaller one
        Cell* baseCell = B0table->getLargestGainCell();
        previousBlock = baseCell->block;
        baseCell->lock();
        lockedCells->insert(baseCell);
        baseCell->switchBlock();
        B0table->remove(baseCell);
        B1table->size+=baseCell->size;
        Net* chosenNet;
        for (int x = 0; x < baseCell->netlist.size(); x++) {
            chosenNet = baseCell->netlist[x];
            chosenNet->updateGains(previousBlock, baseCell, B0table, B1table, cutsetSize);
        }
    }
}

void freeLockedCells(gainTable* B0table, gainTable* B1table, unordered_set<Cell*>* lockedCells) {
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
        Btable->numCells++;
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
    int totalnumCells = (B0table->numCells + B1table->numCells);
    //cout << "total num of cells is " << totalnumCells << "\r\n";
    while (lockedCells->size() < totalnumCells) {//B0 is smaller one
        //cout << "count is " << count << "\r\n";
        nodecnt++;
        count++;
        if (B0table->size >= (totalCells * ratio) && B0table->numCells > 0) {
            previousTable = B0table;
            nextTable = B1table;
        }
        else if (B1table->size > (totalCells * float(1-ratio)) && B1table->numCells > 0) {
            previousTable = B1table;
            nextTable = B0table;
        }
        else {
            if (B0table->numCells == 0 || B1table->numCells == 0)
                cout << "one of the tables has 0 cells in it\r\n";//happens not too often
            break;
        }
        baseCell = previousTable->getLargestGainCell();
        previousBlock = baseCell->block;
        baseCell->lock();
        lockedCells->insert(baseCell);
        baseCell->switchBlock();
        previousTable->remove(baseCell);
        nextTable->size+=baseCell->size;
        //nextTable->numCells++;
        for (auto* x : baseCell->netlist)
            x->updateGains(previousBlock, baseCell, B0table, B1table, cutsetSize);
        if (cutsetSize < bestCutsetSize)
            bestCutsetSize = cutsetSize;
    }
    nodecnt = 0;
}

void freeLockedCellsByCellNum(gainTable* B0table, gainTable* B1table, unordered_set<Cell*>* lockedCells) {
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
        //Btable->numCells++;
    }
    lockedCells->clear();
}

void firstPassByCellNum(gainTable* B0table, gainTable* B1table, unordered_set<Cell*>* lockedCells, float ratio, float totalCells, int& cutsetSize, bool emptyTable) {
    int previousBlock;
    gainTable* previousTable;
    gainTable* nextTable;
    if (emptyTable == 1) {
        previousTable = B0table;
        nextTable = B1table;
    }
    else {
        previousTable = B1table;
        nextTable = B0table;
    }
    while (previousTable->numCells > (totalCells * ratio)) {
        Cell* baseCell = previousTable->getLargestGainCell();
        previousBlock = baseCell->block;
        baseCell->lock();
        lockedCells->insert(baseCell);
        baseCell->switchBlock();
        previousTable->remove(baseCell);
        nextTable->numCells++;
        Net* chosenNet;
        for (int x = 0; x < baseCell->netlist.size(); x++) {
            chosenNet = baseCell->netlist[x];
            chosenNet->updateGains(previousBlock, baseCell, B0table, B1table, cutsetSize);
        }

    }
}

void PassByCellNum(gainTable* B0table, gainTable* B1table, unordered_set<Cell*>* lockedCells, float ratio, int& cutsetSize, int& bestCutsetSize, float totalCells) {
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
        if (previousTable->numCells == 0)
            break;
        baseCell = previousTable->getLargestGainCell();
        previousBlock = baseCell->block;
        baseCell->lock();
        lockedCells->insert(baseCell);
        baseCell->switchBlock();
        previousTable->remove(baseCell);
        nextTable->numCells++;
        for (auto* x : baseCell->netlist)
            x->updateGains(previousBlock, baseCell, B0table, B1table, cutsetSize);
        if (cutsetSize < bestCutsetSize)
            bestCutsetSize = cutsetSize;
    }
}

vector<unordered_set<Cell*>> FMalgoOptimized(float ratio, unordered_set<Cell*>* cellset) {

    
    int count = 0;
    Cell* chosenCell;

    unordered_set<Net*> allNets;
    for (auto x : *cellset) {
        for (auto y : x->netlist) {
            allNets.insert(y);
        }
    }
    unordered_set<Net*> saved_allNets;
    Net* copyNet;

    //UPDATE GAINS FOR ALL CELLS BASED OFF OF BEING IN THE SAME PARTITION, 
    //AND FIND LARGEST AMOUNT OF NETS FOR A CELL, HENCE THE HIGHEST AVAILABLE GAIN
    int largestGain = 0;
    int size = 0;
    for (auto x : *cellset) {
        chosenCell = x;
        chosenCell->initGain();
        size = -1 * chosenCell->gain;
        if (size > largestGain)
            largestGain = size;
    }

    //CREATE HASH TABLES AND INSERT ALL NODES INTO B0
    gainTable* B0table = new gainTable(largestGain);
    gainTable* B1table = new gainTable(largestGain);
    B0table->currentHighestGain = largestGain;
    B1table->currentHighestGain = largestGain;
    for (auto x : *cellset) {
        B0table->insert(x);
    }
    
    float totalCellArea = B0table->size;
    float totalCells = B0table->numCells;
    //cout << "total sieze is " << totalCells << "\r\n";

    //PERFORM FIRST PASS TO MAKE THE FIRST PARTITION
    int cutsetSize = 0;
    unordered_set<Cell*> lockedCells{};
    //cout << "begin first pass\r\n";
    firstPass(B0table, B1table, &lockedCells, ratio, totalCellArea, cutsetSize);
    
    //FREE LOCKED CELLS AND UPDATE POSITIONING FOR ALL LOCKED CELLS
    freeLockedCells(B0table, B1table, &lockedCells);
    bool partByCellnum = false;
    bool emptyTable = 0;
    if (B1table->numCells == 0 && B0table->numCells == 0)
        cout << "\r\n";
    if (B1table->numCells == 0) {
        partByCellnum = true;
        emptyTable = 1;
    }
    if (B0table->numCells == 0) {
        partByCellnum = true;
        emptyTable = 0;
    }
    if (partByCellnum) {
        firstPassByCellNum(B0table, B1table, &lockedCells, ratio, totalCells, cutsetSize, emptyTable);
        freeLockedCellsByCellNum(B0table, B1table, &lockedCells);
    }

    //PERFORM PASSES UNTILL CUTSETSIZE NO LONGER IMPROVES
    count = 0;
    int bestCutesetSize = cutsetSize;
    int previousbestCutsetSize = bestCutesetSize + 1;
    while (bestCutesetSize < previousbestCutsetSize) {
        previousbestCutsetSize = bestCutesetSize;
        
        if (previousbestCutsetSize == 0)
            break;
        if (partByCellnum) {
            PassByCellNum(B0table, B1table, &lockedCells, ratio, cutsetSize, bestCutesetSize, totalCells);
            freeLockedCellsByCellNum(B0table, B1table, &lockedCells);
        }
        else {
            Pass(B0table, B1table, &lockedCells, ratio, cutsetSize, bestCutesetSize, totalCellArea);
            freeLockedCells(B0table, B1table, &lockedCells);
        }
        //cout << "best cutsetsize for this pass is " << bestCutesetSize << "\r\n";
        count++;
    }
    cout << "ratio is " << ratio << " the best cutsetSize is " << previousbestCutsetSize << " number of passes is "
        << count << "\r\n";

    B0size = B0table->size;

    unordered_set<Cell*> B0cells;
    unordered_set<Cell*> B1cells;
    vector<unordered_set<Cell*>> cellLists;

    unordered_set<Net*> B0nets;
    unordered_set<Net*> B1nets;
    for (auto x : allNets) {
        copyNet = new Net;
        copyNet->B0cellList = x->B0cellList;
        B0nets.insert(copyNet);
    }
    for (auto x : allNets) {
        copyNet = new Net;
        copyNet->B0cellList = x->B1cellList;
        B1nets.insert(copyNet);
        delete(x);
    }
    for (auto x : *cellset) {
        x->netlist.clear();
        x->block = 0;
    }
    for (auto x : B0nets) {
        for (auto y : x->B0cellList) {
            y->netlist.push_back(x);
            B0cells.insert(y);
        }
    }
    cout << "cellnum is " << B0cells.size() << "\r\n";
    for (auto x : B1nets) {
        for (auto y : x->B0cellList) {
            y->netlist.push_back(x);
            B1cells.insert(y);
        }
    }
    cout << "cellnum is " << B1cells.size() << "\r\n";
    //saved_allNets.clear();
    cellLists.push_back(B0cells);
    cellLists.push_back(B1cells);
    B0cells.clear();
    B1cells.clear();
    allNets.clear();
    delete(B0table);
    delete(B1table);
    return cellLists;

}

unordered_set<Net*> original_nets;
unordered_set<Cell*> pads;
unordered_set<Cell*> traversedCells;
vector<Cell*> newsort;

unordered_set<Cell*> getCells(string sizefile, string netfile) {

    unordered_map<string, Cell*> cellmap;
    unordered_set<Cell*> cellset;
    string name;
    int cellsize;
    ifstream infile;
    infile.open(sizefile);
    if (!infile)
        cout << "could not open file\r\n";
    while (infile >> name >> cellsize) {
        Cell* newcell = new Cell(0, cellsize, name);
        cellmap.insert(make_pair(name, newcell));
        cellset.insert(newcell);
        /*if (name.substr(0, 1)._Equal("p")) {
            newcell->pad = true;
            pads.insert(newcell);
        }*/
    }
    infile.close();

    //CREATE ALL NETS AND PLUG NET TO CELLS AND CELLS TO NETS
    //cout << "create new nets\r\n";
    int count = 0;
    string startNet;
    Cell* chosenCell;
    Net* currentNet = nullptr;
    Net* currentNet2 = nullptr;

    infile.open(netfile);
    if (!infile)
        cout << "could not open file\r\n";
    std::string line;
    while (std::getline(infile, line))
    {
        if (count >= 5) {
            std::istringstream iss(line);
            if (!(iss >> name >> startNet)) { break; } // error
            //cout << "name is " << name << "\r\n";
            if (startNet._Equal("s")) {
                
                currentNet = new Net();
                currentNet2 = new Net();
                original_nets.insert(currentNet2);
            }
            
            chosenCell = cellmap.at(name);
            chosenCell->addNet(currentNet);
            currentNet->addNode(chosenCell);
            currentNet2->addNode(chosenCell);
        }
        else
            count++;
    }
    infile.close();
    cout << "size of output pads is " << pads.size() << "\r\n";


    return cellset;

}

vector<unordered_set<Cell*>> recur_vert_cuts0(unordered_set<Cell*>* cellset) {
    vector<unordered_set<Cell*>> cellLists;
    vector<unordered_set<Cell*>> lists1;
    vector<unordered_set<Cell*>> lists2;
    unordered_set<Cell*> buffer;
    int cnt = 0;
    int tempsize = 0;
    
    if (cellset->size() <= 2) {
        for (auto x : *cellset) {
            if (cnt == 1)
                x->xpos += (tempsize / hmin);
            buffer.insert(x);
            cellLists.push_back(buffer);
            buffer = {};
            cnt++;
            tempsize = x->size;
        }
    }
    else {
        cout << "before partition\r\n";
        cellLists = FMalgoOptimized(0.5, cellset);
        cout << "after partition\r\n";
        for (auto x : cellLists[1])
            x->xpos += (B0size / hmin);
        
        lists1 = recur_vert_cuts0(&cellLists[0]);
        lists2 = recur_vert_cuts0(&cellLists[1]);

        //cellLists.clear();
        
        cellLists = lists1;
        lists1.clear();
        for (auto x : lists2) 
            cellLists.push_back(x);
        lists2.clear();
        lists1.clear();
    }
    return cellLists;
}

vector<unordered_set<Cell*>> recur_horz_cuts0(unordered_set<Cell*>* cellset, int h) {
    vector<unordered_set<Cell*>> cellLists;
    vector<unordered_set<Cell*>> lists1;
    vector<unordered_set<Cell*>> lists2;
    int ypos = (h/2);
    unordered_set<Net*> the_nets = {};
    if (h < 2 * hmin) {
        cellLists.push_back(*cellset);
    }
    else {
        cout << "before partition\r\n";
        cellLists = FMalgoOptimized(0.5, cellset);
        cout << "after partition\r\n";

        if (cellLists[0].size() > 1)
            lists1 = recur_horz_cuts0(&cellLists[0], ypos);
        else
            lists1.push_back(cellLists[0]);
        if (cellLists[1].size() > 1)
            lists2 = recur_horz_cuts0(&cellLists[1], ypos);
        else
            lists2.push_back(cellLists[1]);

        cellLists.clear();

        cellLists = lists1;
        lists1.clear();
        
        for (auto x : lists2) 
            cellLists.push_back(x);
        lists2.clear();
        
    }
    return cellLists;
}

vector<unordered_set<Cell*>> partitionByPriority(unordered_set<Cell*>* cellset, vector<Cell*> sortedCells, float ratio) {
    vector<unordered_set<Cell*>> cellLists;
    unordered_set<Cell*> set1;
    unordered_set<Cell*> set2;
    //vector<Cell*> newsort;
    int area = 0;
    
    for (int x = 0; x < sortedCells.size(); x++) {
        if (area < float(B0size * ratio)) {
            newsort.push_back(sortedCells[x]);
            set1.insert(sortedCells[x]);
            area += sortedCells[x]->size;
        }
        else {
            set2.insert(sortedCells[x]);
        }
    }
    B0size = area;

    unordered_set<Net*> allNets;
    for (auto x : *cellset) {
        for (auto y : x->netlist) {
            allNets.insert(y);
        }
    }

    unordered_set<Net*> B0nets;
    unordered_set<Net*> B1nets;
    for (auto net : allNets) {
        Net* net1 = new Net;
        Net* net2 = new Net;
        for (auto cell : net->B0cellList) {
            if (set1.count(cell))
                net1->addNode(cell);
            else
                net2->addNode(cell);
        }
        B0nets.insert(net1);
        B1nets.insert(net2);
        delete(net);
    }

    for (auto x : *cellset) {
        x->netlist.clear();
        x->block = 0;
    }

    for (auto x : B0nets) {
        for (auto y : x->B0cellList) {
            y->netlist.push_back(x);
        }
    }
    for (auto x : B1nets) {
        for (auto y : x->B0cellList) {
            y->netlist.push_back(x);
        }
    }
    cellLists.push_back(set1);
    cellLists.push_back(set2);
    return cellLists;
}


vector<unordered_set<Cell*>> recur_horz_cuts_priority0(unordered_set<Cell*>* cellset, int h, vector<Cell*> sortedCells) {
    vector<unordered_set<Cell*>> cellLists;
    vector<unordered_set<Cell*>> lists1;
    vector<unordered_set<Cell*>> lists2;
    int ypos = (h / 2);
    unordered_set<Net*> the_nets = {};
    if (h < 2 * hmin) {
        cellLists.push_back(*cellset);
    }
    else {
        cout << "before partition\r\n";
        cellLists = partitionByPriority(cellset, sortedCells, 0.5);
        cout << "after partition\r\n";

        if (cellLists[0].size() > 1)
            lists1 = recur_horz_cuts_priority0(&cellLists[0], ypos, newsort);
        else
            lists1.push_back(cellLists[0]);
        if (cellLists[1].size() > 1)
            lists2 = recur_horz_cuts0(&cellLists[1], ypos);
        else
            lists2.push_back(cellLists[1]);


        cellLists = lists1;
        lists1.clear();
        newsort.clear();

        for (auto x : lists2)
            cellLists.push_back(x);
        lists2.clear();

    }
    return cellLists;
}

vector<unordered_set<Cell*>> recur_vert_cuts_priority0(unordered_set<Cell*>* cellset, vector<Cell*> sortedCells) {
    vector<unordered_set<Cell*>> cellLists;
    vector<unordered_set<Cell*>> lists1;
    vector<unordered_set<Cell*>> lists2;
    
    unordered_set<Cell*> buffer;
    int cnt = 0;
    int tempsize = 0;

    if (cellset->size() <= 2) {
        for (auto x : *cellset) {
            if (cnt == 1)
                x->xpos += (tempsize / hmin);
            buffer.insert(x);
            cellLists.push_back(buffer);
            buffer = {};
            cnt++;
            tempsize = x->size;
        }
    }
    else {
        cout << "before partition\r\n";
        cellLists = partitionByPriority(cellset, sortedCells, 0.5);
        cout << "after partition\r\n";
        for (auto x : cellLists[1])
            x->xpos += (B0size / hmin);

        if (cellLists[0].size() > 1)
            lists1 = recur_vert_cuts_priority0(&cellLists[0], newsort);
        else
            lists1.push_back(cellLists[0]);
        if (cellLists[1].size() > 1)
            lists2 = recur_vert_cuts0(&cellLists[1]);
        else
            lists2.push_back(cellLists[1]);

        cellLists = lists1;
        lists1.clear();
        for (auto x : lists2)
            cellLists.push_back(x);
        lists2.clear();
    }
    return cellLists;
}

vector<unordered_set<Cell*>> recur_sq_cuts(unordered_set<Cell*>* cellset, bool vert) {
    vector<unordered_set<Cell*>> cellLists;
    vector<unordered_set<Cell*>> lists1;
    vector<unordered_set<Cell*>> lists2;
    unordered_set<Cell*> buffer;
    int cnt = 0;
    int tempsize = 0;

    if (cellset->size() <= 2) {
        for (auto x : *cellset) {
            if (cnt == 1) {
                if (vert)
                    x->xpos += int(sqrt(tempsize));
                else
                    x->ypos += int(sqrt(tempsize));
            }
            buffer.insert(x);
            cellLists.push_back(buffer);
            buffer = {};
            cnt++;
            tempsize = x->size;
        }
    }
    else {
        cout << "before partition\r\n";
        cellLists = FMalgoOptimized(0.5, cellset);
        cout << "after partition\r\n";
        for (auto x : cellLists[1]) {
            if (vert)
                x->xpos += int(sqrt(B0size));
            else
                x->ypos += int(sqrt(B0size));
        }

        lists1 = recur_sq_cuts(&cellLists[0], !vert);
        lists2 = recur_sq_cuts(&cellLists[1], !vert);

        cellLists.clear();

        cellLists = lists1;
        for (auto x : lists2) {
            //for (auto y : x)
                //y->xpos += xpos;
            cellLists.push_back(x);
        }
        lists2.clear();
        lists1.clear();
    }
    return cellLists;
}

vector<unordered_set<Cell*>> recur_sq_cuts_priority(unordered_set<Cell*>* cellset, bool vert, vector<Cell*> sortedCells) {
    vector<unordered_set<Cell*>> cellLists;
    vector<unordered_set<Cell*>> lists1;
    vector<unordered_set<Cell*>> lists2;
    unordered_set<Cell*> set1;
    unordered_set<Cell*> set2;
    vector<Cell*> newsort;
    unordered_set<Cell*> buffer;
    int cnt = 0;
    int tempsize = 0;
    int area = 0;
    int numcells = 0;
    if (cellset->size() <= 2) {
        for (auto x : *cellset) {
            if (cnt == 1) {
                if (vert)
                    x->xpos += int(sqrt(tempsize));
                else
                    x->ypos += int(sqrt(tempsize));
            }
            buffer.insert(x);
            cellLists.push_back(buffer);
            buffer = {};
            cnt++;
            tempsize = x->size;
        }
    }
    else {
        cout << "before partition\r\n";
        cellLists = partitionByPriority(cellset, sortedCells, 0.5);
        cout << "after partition\r\n";
        for (auto x : cellLists[1]) {
            if (vert)
                x->xpos += int(sqrt(B0size));
            else
                x->ypos += int(sqrt(B0size));
        }

        lists1 = recur_sq_cuts_priority(&cellLists[0], !vert, newsort);
        lists2 = recur_sq_cuts(&cellLists[1], !vert);

        cellLists = lists1;
        for (auto x : lists2) {
            //for (auto y : x)
                //y->xpos += xpos;
            cellLists.push_back(x);
        }
        lists2.clear();
        lists1.clear();
    }
    return cellLists;
}

int vertcutcnt = 0;
int width = 0;
vector<unordered_set<Cell*>> recur_horz_cuts_priority1(unordered_set<Cell*>* cellset) {
    cout << "recur is " << vertcutcnt << " deep" << "\r\n";
    vertcutcnt++;
    vector<unordered_set<Cell*>> cellLists;
    vector<unordered_set<Cell*>> lists1;
    float h = float(B0size / width);
    float ratio = (width * (h - hmin))/B0size;
    cout << "current height is " << h << " num of cells is " << cellset->size() << "\r\n";
    unordered_set<Net*> the_nets = {};
    if (h < 2 * hmin) {
        //cellLists = new vector<unordered_set<Cell*>>;
        cellLists.push_back(*cellset);
    }
    else {
        cout << "before partition\r\n";
        cellLists = FMalgoOptimized(ratio, cellset);
        cout << "after partition\r\n";

        cout << "cellLists size is " << cellLists.size() << "\r\n";
        if (cellLists[0].size() > 1)
            lists1 = recur_horz_cuts_priority1(&cellLists[0]);
        else
            lists1.push_back(cellLists[0]);
        

        cellLists.erase(cellLists.begin());//removes [0]


        lists1.push_back(cellLists[0]);
        cellLists = lists1;
        lists1.clear();

    }
    cout << "recur is " << vertcutcnt << " deep" << "\r\n";
    vertcutcnt--;
    return cellLists;
}

vector<unordered_set<Cell*>> recur_vert_cuts_priority1(unordered_set<Cell*>* cellset) {
    vector<unordered_set<Cell*>> cellLists;
    vector<unordered_set<Cell*>> lists1;
    unordered_set<Cell*> buffer{};
    Cell* chosenCell;
    int cnt = 0;
    int tempsize = 0;

    int largestGain = 0;
    int gain = 0;
    int pos = 0;
    unordered_map<int, unordered_set<Cell*>> sortmap;
    for (auto cell : *cellset) {
        chosenCell = cell;
        chosenCell->initGain();
        gain = -1 * cell->gain;
        if (!sortmap.count(gain))
            sortmap[gain] = buffer; 
        sortmap[gain].insert(cell);
        if (gain > largestGain) 
            largestGain = gain;
        //cout << "gain is " << cell->gain << "\r\n";
        //cout << "bucket count is " << sortmap.size() << "\r\n";
    }
    for (int x = 0; x <= largestGain; x++) {
        if (sortmap.count(x)) {
            //cout << "length of bucket is " << sortmap[x].size() << "\r\n";
            for (auto cell : sortmap[x]) {
                cell->xpos = pos;
                pos += (cell->size) / hmin;
                buffer.insert(cell);
                cellLists.push_back(buffer);
                buffer.clear();
            }
        }
    }
    return cellLists;
        
    
}

int wirelengthX(Net* n, int width) {
    int minX = width;
    int maxX = 0;
    for (auto x : n->B0cellList) {
        if (x->xpos < minX)
            minX = x->xpos;
        if (x->xpos > maxX)
            maxX = x->xpos;
    }
    return (maxX - minX);
}

int wirelengthY(Net* n, int height) {
    int minY = height;
    int maxY = 0;
    for (auto x : n->B0cellList) {
        if (x->ypos < minY)
            minY = x->ypos;
        if (x->ypos > maxY)
            maxY = x->ypos;
    }
    return (maxY - minY);
}

vector<Cell*> sortByPriority(unordered_set<Cell*>* cellset) {
    unordered_map<int, unordered_set<Cell*>> sortmapP;
    vector<Cell*> sortedCells;
    Cell* chosenCell;
    unordered_set<Cell*> buffer{};
    int priority = 0;
    int highestPriority = 0;
    for (auto cell : *cellset) {
        chosenCell = cell;
        chosenCell->initGain();
        priority = cell->size;
        if (!sortmapP.count(priority))
            sortmapP[priority] = buffer;
        sortmapP[priority].insert(cell);
        if (priority > highestPriority)
            highestPriority = priority;
    }
    for (int x = highestPriority; x >= 0; x--) {
        if (sortmapP.count(x)) {
            for (auto cell : sortmapP[x]) {
                sortedCells.push_back(cell);
            }
        }
    }
    return sortedCells;
}

int main() {

    /*  Testing Directions:
            -use horzcuts0 and vertcuts0 for bisection 
            -sqcuts for quadtrature
            -horzcuts1 and vertcuts0 for slice bi section
            -horzcuts1 and vertcuts1 for cut oriented
    */


    auto start = high_resolution_clock::now();
    srand(time(NULL));

    float ratio[] = { 0.5, 0.45, 0.4, 0.35 };
    string sizefile[] = { "C:\\temp\\ibm01.are", "C:\\temp\\ibm02.are", "C:\\temp\\ibm03.are", "C:\\temp\\ibm04.are", "C:\\temp\\ibm05.are", "C:\\temp\\ibm06.are"
                        , "C:\\temp\\ibm07.are" , "C:\\temp\\ibm08.are" , "C:\\temp\\ibm09.are" , "C:\\temp\\ibm10.are" , "C:\\temp\\ibm11.are" , "C:\\temp\\ibm12.are"
                        , "C:\\temp\\ibm13.are" , "C:\\temp\\ibm14.are" , "C:\\temp\\ibm15.are" , "C:\\temp\\ibm16.are" , "C:\\temp\\ibm17.are" , "C:\\temp\\ibm18.are" };
    string netfile[] = { "C:\\temp\\ibm01.net", "C:\\temp\\ibm02.net", "C:\\temp\\ibm03.net", "C:\\temp\\ibm04.net", "C:\\temp\\ibm05.net", "C:\\temp\\ibm06.net"
                        , "C:\\temp\\ibm07.net" , "C:\\temp\\ibm08.net" , "C:\\temp\\ibm09.net" , "C:\\temp\\ibm10.net" , "C:\\temp\\ibm11.net" , "C:\\temp\\ibm12.net"
                        , "C:\\temp\\ibm13.net" , "C:\\temp\\ibm14.net" , "C:\\temp\\ibm15.net" , "C:\\temp\\ibm16.net" , "C:\\temp\\ibm17.net" , "C:\\temp\\ibm18.net" };

    unordered_set<Cell*> the_set;
    the_set = getCells(sizefile[0], netfile[0]);
    float totalArea = 0;
    for (auto x : the_set) {
        totalArea += x->size;
    }
    B0size = totalArea;
    int h = int(sqrt(totalArea));
    width = h;
    cout << "total height is " << h << "\r\n";

    int minsize = 100000;
    for (auto cell : the_set) {
        if (cell->size < minsize)
            if(cell->size!=0)
                minsize = cell->size;
    }
    cout << "min cellsize is " << minsize << "\r\n";

    vector<Cell*> sortedCells;
    //sortedCells = sortByPriority(&the_set);

    //vector<unordered_set<Cell*>> result = recur_sq_cuts_priority(&the_set, true, sortedCells);
    vector<unordered_set<Cell*>> result = recur_horz_cuts0(&the_set, h);
    //vector<unordered_set<Cell*>> result = recur_horz_cuts_priority1(&the_set);
    for (int x = 0; x < result.size(); x++) {             //COMMENT THIS FOR-LOOP FOR QUADRATURE
        for (auto cell : result[x])
            cell->ypos = (hmin * x);
    }
    cout << "length of list of cell sets (aka partitions) is " << result.size() << "\r\n";

    unordered_set<int> sizes;

    cout << "*******************************************************************************************" << "\r\n";
    vector<unordered_set<Cell*>> result2;
    vector<unordered_set<Cell*>> temp;
    for (auto x : result) {                               //COMMENT THIS FOR-LOOP FOR QUADRATURE PLACEMENT
        //sortedCells = sortByPriority(&x);
        temp = recur_vert_cuts0(&x);
        //temp = recur_vert_cuts_priority1(&x);
        for (auto y : temp)
            result2.push_back(y);
    }
    temp.clear();
    cout << "length of list of cell sets (aka partitions) is " << result2.size() << "\r\n";
    
    
    
    int totalWireLengthX = 0;
    int totalWireLengthY = 0;
    for (auto x : original_nets) {
        totalWireLengthX += wirelengthX(x, h);
        totalWireLengthY += wirelengthY(x, h);
    }
    cout << "toal x length is " << totalWireLengthX << "\r\n";
    cout << "toal y length is " << totalWireLengthY << "\r\n";

    int minwidth = minsize / hmin;
    int count = 0;
    int maxrowlength = 0;
    int maxcollength = result.size();
    unordered_set<int> positions;
    for (auto row : result) {
        int previouspos = 0;
        int previouswidth = 0;
        if (row.size() > maxrowlength)
            maxrowlength = row.size();
        for (auto cell : row) {
            positions.insert(cell->xpos);
            cell->xpos = ceil(float((cell->xpos) / minwidth));
            cell->ypos = (cell->ypos) / hmin;
            //positions.insert(cell->xpos);
            //cout << "xpos is " << cell->xpos << "\r\n";
        }
        if (positions.size() != row.size()) {
            count++;
            cout << "error found \r\n";
        }
        positions.clear();
        //count++;
    }
    cout << "errors found is " << count << "\r\n";
    cout << "rows is " << maxcollength << "\r\n";
    cout << "cols is " << maxrowlength << "\r\n";


    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << " execution time(ms) is " << duration.count() << "\r\n";
    cout << "done";
    return 0;
}
