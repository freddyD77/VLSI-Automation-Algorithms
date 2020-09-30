
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

int hmin = 16;

class Net;
class netList;
class Cell;

class Vertex
{
    // Access specifier 
public:

    // Data Members 
    int positionx;
    int positiony;
    vector<int> neighbors;
    Cell* cell;
    bool blocked;
    int level=-1;
    int nets;

    Vertex(int posx, int posy) {
        positionx = posx;
        positiony = posy;
        nets = 0;
        blocked = false;
    }

    vector<int> getPos() {
        vector<int> pos;
        pos.push_back(positionx);
        pos.push_back(positiony);
        return pos;
    }

};

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

};

class Grid
{
    // Access specifier 
public:

    // Data Members 
    int rows;
    int cols;
    int currentfeeds;
    vector<vector<Vertex*>> vertices;
    unordered_set<Vertex*> feedthroughs;
    unordered_set<Vertex*> traversed;

    Grid(int numRows, int numCols, unordered_set<Cell*>* cellset) {
        cols = numCols;
        rows = numRows;
        vector<Vertex*> verticerow;
        for (int y = 0; y < numRows; y++) {
            verticerow.clear();
            for (int x = 0; x < numCols; x++) {
                Vertex* vertice = new Vertex(x, y);
                verticerow.push_back(vertice);
            }
            vertices.push_back(verticerow);
        }

        for (auto* cell : *cellset) {
            vertices[cell->xpos][cell->ypos]->blocked = true;//block every vertice that has a cell
            vertices[cell->xpos][cell->ypos]->cell = cell;//put cell in vertex
        }

        currentfeeds = 0;
    }

    void insertFeedThrough(int feeds) {//adds feed throughs symmetrically for each row
        int feedinc = cols / feeds;
        int feedpos = cols / (2 * feeds);
        for (auto row : vertices) {
            while (feedpos < cols) {
                row[feedpos]->blocked = false;//unblock to make a feed through
                feedthroughs.insert(row[feedpos]);
                feedpos += feedinc;
            }
        }
    }

    vector<Vertex*> getNeighbors(vector<int> coord) {
        vector<Vertex*> verts;
        int up = coord[1] + 1;
        int down = coord[1] - 1;
        int right = coord[0] + 1;
        int left = coord[0] - 1;
        if (!up > rows)
            verts.push_back(vertices[coord[0]][up]);
        if (!down < 0)
            verts.push_back(vertices[coord[0]][down]);
        if (!right > cols)
            verts.push_back(vertices[right][coord[1]]);
        if (!left < 0)
            verts.push_back(vertices[left][coord[1]]);
        return verts;
    }

    Vertex* getNeighborInDir(Vertex* vert, vector<int> coord);

};

class Net {
public:
    vector<Cell*> cellList;
    unordered_set<Net*> netsBelow;
    Net() {}



    void addNode(Cell* cellin) {
        cellList.push_back(cellin);

    }

    void printB0list() {
        for (int x = 0; x < cellList.size(); x++) {
            cout << cellList.at(x) << " ";
        }
    }


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

unordered_set<Net*> original_nets;

unordered_set<Cell*> getCells(string sizefile, string netfile, string posfile) {

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
    }
    infile.close();

    //CREATE ALL NETS AND PLUG NET TO CELLS AND CELLS TO NETS
    //cout << "create new nets\r\n";
    int count = 0;
    string startNet;
    Cell* chosenCell;
    Net* currentNet = nullptr;

    infile.open(netfile);
    if (!infile)
        cout << "could not open file\r\n";
    std::string line;
    while (std::getline(infile, line))
    {
        if (count >= 5) {
            std::istringstream iss(line);
            if (!(iss >> name >> startNet)) { break; } // error
            if (startNet._Equal("s")) {

                currentNet = new Net();
                original_nets.insert(currentNet);
            }

            chosenCell = cellmap.at(name);
            chosenCell->addNet(currentNet);
            currentNet->addNode(chosenCell);
        }
        else
            count++;
    }
    infile.close();

    string xpos;
    string ypos;
    infile.open(posfile);//gets the position data from placement results
    if (!infile)
        cout << "could not open file\r\n";
    std::string line2;
    while (std::getline(infile, line2))
    {
        std::istringstream iss(line2);
        if (!(iss >> name >> xpos >> ypos)) { break; } // error
        
        chosenCell = cellmap.at(name);
        chosenCell->xpos = stoi(xpos);
        chosenCell->ypos = stoi(ypos);
    }
    infile.close();


    return cellset;

}

int getDirectionX(vector<int> coord1, vector<int> coord2) {
    if (coord2[0] > coord1[0])
        return 2;
    if (coord2[0] == coord1[0])
        return 1;
    if (coord2[0] < coord1[0])
        return 0;
}

int getDirectionY(vector<int> coord1, vector<int> coord2) {
    if (coord2[1] > coord1[1])
        return 2;
    if (coord2[1] == coord1[1])
        return 1;
    if (coord2[1] < coord1[1])
        return 0;
}

Vertex* Grid::getNeighborInDir(Vertex* vert, vector<int> Endcoords) {
    vector<int> coord = vert->getPos();
    vector<Vertex*> neighbors = getNeighbors(coord);
    for (auto neighb : neighbors) {
        if (getDirectionX(coord, Endcoords) == getDirectionX(neighb->getPos(), Endcoords) ||
            getDirectionY(coord, Endcoords) == getDirectionY(neighb->getPos(), Endcoords)) {
            return neighb;
        }
    }

    return nullptr;
}

void Retrace(int startx, int starty, int endx, int endy, Grid* grid) {
    vector<Vertex*> neighbors;
    Vertex* currentVert = grid->vertices[endx][endy];
    Vertex* startVert = grid->vertices[startx][starty];
    startVert->level = 0;
    Vertex* nextVert=nullptr;
    vector<int> coord;
    int currentLevel;
    while (currentVert != startVert) {
        coord.push_back(currentVert->positionx);
        coord.push_back(currentVert->positiony);
        neighbors = grid->getNeighbors(coord);
        currentLevel = currentVert->level;
        for (auto neighbor : neighbors) {
            if (neighbor->level < currentLevel && neighbor->level != -1) {//-1 is levels of non travrsed vertices
                currentLevel = neighbor->level;
                nextVert = neighbor;
            }
        }//result is nextVert is the next vertice (with lowest level) to retrace
        currentVert = nextVert;
        grid->traversed.insert(currentVert);//add traversed vertices to set
    }
}

bool Soukup(int startx, int starty, int endx, int endy, Grid* grid) {
    vector<vector<int>> previousVertsXandY;//contains (x,y) coordinates of vertices
    vector<vector<int>> newVertsXandY;
    vector<Vertex*> neighbors;
    vector<int> coords;
    vector<int> Endcoords;
    coords.push_back(startx);
    coords.push_back(starty);
    Endcoords.push_back(endx);
    Endcoords.push_back(endy);
    previousVertsXandY.push_back(coords);
    int temp = 1;
    int *xcoord;
    int *ycoord;
    bool path_exists = false;
    while (previousVertsXandY.size() != 0){
        for (auto coord : previousVertsXandY) {
            neighbors = grid->getNeighbors(coord);
            for (auto neighbor : neighbors) {
                xcoord = &neighbor->positionx;
                ycoord = &neighbor->positiony;
                if (*xcoord == endx && *ycoord == endy) {
                    grid->vertices[endx][endy]->level = temp;
                    path_exists = true;
                    return path_exists;
                }
                if (!neighbor->blocked) {
                    if (getDirectionX(coord, Endcoords) == getDirectionX(neighbor->getPos(), Endcoords) ||
                        getDirectionY(coord, Endcoords) == getDirectionY(neighbor->getPos(), Endcoords)) {
                        neighbor->level = temp;
                        temp++;
                        previousVertsXandY.push_back(neighbor->getPos());
                        while (!grid->getNeighborInDir(neighbor, Endcoords)->blocked) {
                            neighbor = grid->getNeighborInDir(neighbor, Endcoords);
                            neighbor->level = temp;
                            temp++;
                            previousVertsXandY.push_back(neighbor->getPos());
                        }
                    }
                    else {
                        neighbor->level = temp;
                        temp++;
                        newVertsXandY.push_back(neighbor->getPos());
                    }
                }
            }
        }
        previousVertsXandY = newVertsXandY;
        newVertsXandY.clear();
    }
    return path_exists;
}



int main() {

    auto start = high_resolution_clock::now();
    srand(time(NULL));

    float ratio[] = { 0.5, 0.45, 0.4, 0.35 };
    string sizefile[] = { "C:\\temp\\ibm01.are", "C:\\temp\\ibm02.are", "C:\\temp\\ibm03.are", "C:\\temp\\ibm04.are", "C:\\temp\\ibm05.are", "C:\\temp\\ibm06.are"
                        , "C:\\temp\\ibm07.are" , "C:\\temp\\ibm08.are" , "C:\\temp\\ibm09.are" , "C:\\temp\\ibm10.are" , "C:\\temp\\ibm11.are" , "C:\\temp\\ibm12.are"
                        , "C:\\temp\\ibm13.are" , "C:\\temp\\ibm14.are" , "C:\\temp\\ibm15.are" , "C:\\temp\\ibm16.are" , "C:\\temp\\ibm17.are" , "C:\\temp\\ibm18.are" };
    string netfile[] = { "C:\\temp\\ibm01.net", "C:\\temp\\ibm02.net", "C:\\temp\\ibm03.net", "C:\\temp\\ibm04.net", "C:\\temp\\ibm05.net", "C:\\temp\\ibm06.net"
                        , "C:\\temp\\ibm07.net" , "C:\\temp\\ibm08.net" , "C:\\temp\\ibm09.net" , "C:\\temp\\ibm10.net" , "C:\\temp\\ibm11.net" , "C:\\temp\\ibm12.net"
                        , "C:\\temp\\ibm13.net" , "C:\\temp\\ibm14.net" , "C:\\temp\\ibm15.net" , "C:\\temp\\ibm16.net" , "C:\\temp\\ibm17.net" , "C:\\temp\\ibm18.net" };
    string posfile[] = { "C:\\temp\\ibm01pos.txt", "C:\\temp\\ibm02pos.txt", "C:\\temp\\ibm03pos.txt", "C:\\temp\\ibm04pos.txt", "C:\\temp\\ibm05pos.txt", "C:\\temp\\ibm06pos.txt"
                        , "C:\\temp\\ibm07pos.txt" , "C:\\temp\\ibm08pos.txt" , "C:\\temp\\ibm09pos.txt" , "C:\\temp\\ibm10pos.txt" , "C:\\temp\\ibm11pos.txt" , "C:\\temp\\ibm12pos.txt"
                        , "C:\\temp\\ibm13pos.txt" , "C:\\temp\\ibm14pos.txt" , "C:\\temp\\ibm15pos.txt" , "C:\\temp\\ibm16pos.txt" , "C:\\temp\\ibm17pos.txt" , "C:\\temp\\ibm18pos.txt" };

    int netlimit = 5;//hyper parameter, number of nets per channel

    unordered_set<Cell*> cellset;
    cellset = getCells(sizefile[0], netfile[0], posfile[0]);
    int minsize = 1000000;
    for (auto cell : cellset) {
        if (cell->size > minsize)
            minsize = cell->size;
    }

    int minwidth = minsize / hmin;
    int count = 0;
    int rows = 0;
    int cols = 0;
    //unordered_set<int> positions;
    for (auto cell : cellset) {
        if (cell->xpos > cols)
            cols = cell->xpos;
        if (cell->ypos > rows)
            cols = cell->ypos;
        cell->xpos = ceil(cell->xpos / minwidth);//formatting posiiton vbase doff grid
        cell->ypos = 2*(cell->ypos / hmin);//formatiing posiution based off grid
    }
    rows = rows / hmin;//number of cell rows
    cols = ceil(cols / minwidth);//number of cell columns
    rows = 2 * rows;//we want 1 row spacing in the grid for soukups algo

    int feeds = 1;
    Grid *grid = new Grid(rows, cols, &cellset);
    grid->insertFeedThrough(feeds);//put atleast 1 feed through at each row

    int startx = 0;
    int starty = 0;
    int endx = 0;
    int endy = 0;
    bool pass;
    for (auto net : original_nets) {
        for (int x = 1; x < net->cellList.size(); x++) {//routes every connection in net
            startx = net->cellList[0]->xpos;
            starty = net->cellList[0]->ypos;
            endx = net->cellList[x]->xpos;
            endy = net->cellList[x]->ypos;
            pass = Soukup(startx, starty, endx, endy, grid);
            if (pass) {
                Retrace(startx, starty, endx, endy, grid);//adds vertices to traversed set
            }
            else {
                feeds *= 2;
                grid->insertFeedThrough(feeds);//double the feed throughs
                x--;//re-do this connection since more feedthroughs habve been added
            }
        }
        for (auto vertex : grid->traversed) {//increment net counts for each traversed vertex
            vertex->nets++;
        }
        grid->traversed.clear();//clear traversed vertices
        for (auto vertex : grid->feedthroughs) {
            if (vertex->nets >= netlimit)
                vertex->blocked = true;//block feedthroughs that have hit net limit
        }

    }
    //done with global routing

    //begin making VCGs for every channel(intermeddiate rows)
    vector<unordered_set<Net*>> VCGs;//list of VCGs
    unordered_set<Net*> VCG;//will hold all nets in the same VCG
    vector<Net*> pair;
    pair.push_back(nullptr);
    pair.push_back(nullptr);
    vector<vector<Net*>> pairs;
    int netIndex = 0;
    bool noMoreNets = false;
    for (int y = 1; y < rows-1; y += 2) {//does the rows that do not have standard cells
        while (1) {//only break when there are no more nets to consider
            noMoreNets = true;
            for (int x = 0; x < cols; x++) {
                if (grid->vertices[x][y+2]->cell->netlist.size() > netIndex) {//gets top net
                    pair[0] = grid->vertices[x][y+2]->cell->netlist[netIndex];
                    VCG.insert(pair[0]);
                    noMoreNets = false;
                }
                if (grid->vertices[x][y]->cell->netlist.size() > netIndex) {//gets bottom net
                    pair[1] = grid->vertices[x][y]->cell->netlist[netIndex];
                    VCG.insert(pair[1]);
                    noMoreNets = false;
                }
                pairs.push_back(pair);
            }
            netIndex++;
            if (noMoreNets)
                break;//no more nets at this channel
            //else
            for (auto netpair : pairs) {//for every pair, add the bottom net to the top' list of bottom nets
                netpair[0]->netsBelow.insert(netpair[1]);//Net objects have the VCG connections themselves
            }
        }
        VCGs.push_back(VCG);//saves the VCG of this row/channel
    }

    //begin making zones
    //project ends here, zone creating and merging were too difficult to implement (and conceptualize)
    //under time constarints


    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << " execution time(ms) is " << duration.count() << "\r\n";
    cout << "done with global routing";
    return 0;
}