/*******************************************************************************
//  1TSPbruteForce.cpp              Author: Ian Nobile
//  Section: 50                     Due Date: 30 August 2021
//  
//  This memory-leak- and bug-free program uses a brute force method to find the
//  minimum cost solution to a Traveling Salesperson Problem by creating an 
//  undirected graph from a .TSP file, generating all permutations of a vector 
//  of destinations, tracing a Hamiltonian path and reporting the total distance
*******************************************************************************/

#include <chrono>   // time the speed of the program
#include <iostream> // print to console
#include <fstream>  // read files
#include <string>   // use as a header buffer for seeking inside a file
#include <vector>   // easy arraying
#include <cmath>    // distance formula
using namespace std;
using namespace std::chrono;


// class declarations
class Node {
public:
    int num;
    float x;
    float y;
    bool visited;   // need this later?

    // constructor:
    Node() {
        num = 0;
        x = 0.0;
        y = 0.0;
        visited = false;
    }
};

class Graph {
public:
    vector<Node> nodes;
};

// function prototypes:
Graph buildGraph(int*, char*);
float distCheck(float, float, float, float);
vector<vector<float>> createLookup(int, const Graph&);
int factorial(int);
void permute(int dimension, vector<int> &destination);


//------------------------------------------------------------------------------
//  Main Function
//------------------------------------------------------------------------------
int main(int argc, char* argv[]){
    auto start = high_resolution_clock::now();  // start timer

    // check if path was passed as arg:
    if (argc == 1) {
        cout << "Please pass the path to the .TSP file as a command line argument" << endl;
        return 1;
    }
    // begin with a friendly greeting:
    cout << "Hello and welcome to the (brute force) Travelling Salesperson Problem Solver" << endl;
    int dimension;
    Graph graph = buildGraph(&dimension, argv[1]);
    if (dimension < 1) { // error handling if .TSP could not be opened
        cout << "Error: file could not be opened. Please check the status of your .TSP file" << endl;
        return -1;
    } else if (dimension == 1) { // if there's only one coordinate inside
        cout << "There is only one node and thus, only one possible route of 0 miles." << endl << endl;
        return 0;
    } else {
        vector<vector<float>> distref = createLookup(dimension, graph);
        
        // build destinations vector
        vector<int> destinations(dimension, 0);
        for (int i = 0;i < dimension;i++) {
            destinations[i] = i;
        }

        int dimFact = factorial(dimension);
        cout << "There are " << dimFact << " possible routes. Better start checking them now..." << endl;
        float  minDist = FLT_MAX;
        cout << "The current shortest distance is " << minDist << ". Surely we can do better:" << endl;

        float newDist = 0;
        int first;
        int next;

        for(int i=0;i<= dimFact;i++){
            for (int col = 0;col < dimension - 1;col++) {
                first = destinations[col];
                next = destinations[col+1];
                newDist += distref[first][next];
                if (newDist > minDist) { break; }
            }
            if (newDist > minDist) { 
                newDist = 0;
                if (i == dimFact - 1) { break; } // lexicographic end -- no need to permute again
                permute(dimension, destinations);
                continue; 
            }
            first = destinations[dimension - 1];
            next = destinations[0];
            newDist += distref[first][next];

            if (newDist < minDist) {
                minDist = newDist;
                cout << endl << "Found a shorter distance!" << endl << "The journey from ";
                for (int col = 0;col < dimension;col++) {
                    cout << destinations[col] + 1;
                    cout << " -> ";
                }
                cout << destinations[0] + 1 << " is only " << minDist << " miles in total!" << endl;
            }
            newDist = 0;
            if (i == dimFact - 1) {break;} // lexicographic end -- no need to permute again
            permute(dimension, destinations);
        }
        cout << "And that's the best that can be done." << endl;

        auto stop = high_resolution_clock::now();   // stop timer
        auto duration = duration_cast<microseconds>(stop - start);  // calculate elapsed time
        cout << endl << "Program execution took " << duration.count() / 1000000.0 << "s in total" << endl; // print speed to console

        return 0;
    }
}


//  function definitions:
//------------------------------------------------------------------------------
//  Read .TSP file, creates nodes from coordinates and combines all in an
//  undirected graph object
//------------------------------------------------------------------------------
Graph buildGraph(int* dimension, char* argv) {
    // open .TSP file in read-only mode:
    ifstream tspfile;
    tspfile.open(argv, ios::in);
    // ensure file exists:
    if (!tspfile.is_open()) {
        Graph graph;
        *dimension = 0;
        return graph;
    }
    // advance buffer to the dimension section:
    string name = "";
    while (name.compare("DIMENSION:") != 0) {
        tspfile >> name;
    }
    tspfile >> *dimension;
    // advance buffer to the coordinates section:
    while (name.compare("NODE_COORD_SECTION") != 0) {
        tspfile >> name;
    }
    // create nodes and push to graph vector
    Graph graph;
    Node newNode = Node();
    for (int i = 0;i < *dimension;i++) {
        tspfile >> newNode.num;
        tspfile >> newNode.x;
        tspfile >> newNode.y;
        //newNode.visited = false;  // for later?
        graph.nodes.push_back(newNode);
    }
    // The graph is now created, and we are finished with the .TSP file
    tspfile.close();
    return graph;
}

//------------------------------------------------------------------------------
//  Returns the distance between two graph nodes using the pythagoran formula 
//  dist = sqrt((x2 - x1)^2 + (y2 - y1)^2)(used when constructing for-
//  reference distance matrix)
//------------------------------------------------------------------------------
float distCheck(float x1, float y1, float x2, float y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

//------------------------------------------------------------------------------
//  Creates distance lookup table by calculating distance between each coodinate
//  and every other
//  (Needs optimising -- some repeated work being done here)
//------------------------------------------------------------------------------
vector<vector<float>> createLookup(int dimension, const Graph &graph) {
    vector<vector<float>> distref(dimension, vector<float>(dimension, FLT_MAX));
    for (int row = 0;row < dimension;row++) {
        for (int col = 0;col < dimension;col++) {
            distref[row][col] = distCheck(graph.nodes[row].x, graph.nodes[row].y, graph.nodes[col].x, graph.nodes[col].y);
        }
    }
    return distref;
}

//------------------------------------------------------------------------------
//  Calculates the factorial of a given number (recursive)
//------------------------------------------------------------------------------
int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

//------------------------------------------------------------------------------
//  Generates the next lexicographical permutation of a vector of destinations;
//  based on Algorithm L by the great Donald Knuth:
//------------------------------------------------------------------------------
void permute(int dimension, vector<int> &destinations) {
    int j = dimension - 2;
    int compare = dimension - 1;

    // find next decreasible element
    while (destinations[j] >= destinations[compare]) {
        j--;
        compare--;
        //if (j == -1) { break; }
    }

    int l = dimension - 1;
    while (destinations[j] >= destinations[l]) {
        l--;
    }
    //swap(j,l):
    int temp = destinations[l];
    destinations[l] = destinations[j];
    destinations[j] = temp;

    //reverse everything from [j+1] on:
    j++;
    l = dimension - 1;
    while (j < l) {
        temp = destinations[l];
        destinations[l] = destinations[j];
        destinations[j] = temp;
        j++;
        l--;
    }
}

