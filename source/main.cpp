#include <iostream>
#include <chrono>
#include "ballstar.h"

using namespace std;


int main(){

    BallStar ball("../source/songs_final.csv");
    // BallStar ball("../source/songs_final500.csv");

    cout << "indexando..." << endl;
    
    
    auto start = chrono::steady_clock::now();
    // ball.buildtree;
    ball.buildtree();
    auto end = chrono::steady_clock::now();

    cout << "La construccion duro:" << endl;
    cout << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms" << endl;
    cout << "La construccion duro:" << endl;
        cout << chrono::duration_cast<chrono::seconds>(end - start).count()
        << " sec" << endl;

    Eigen::MatrixXd nn;
    int k;
    // cin>> k;

    cout<<ball.data.row(4)<<endl;

    cout<< "K vecinos (digita K): ";
    cin >> k;
    while (k > 0)
    {
        start = chrono::steady_clock::now();
        ball.searchKNN( nn, ball.data.row(4), k);
        end = chrono::steady_clock::now();
        cout<<"vecinos: " <<endl<<nn<<endl;
        cout << "La busqueda de los vecinos mas cercanos duro:" << endl;
        cout << chrono::duration_cast<chrono::nanoseconds>(end - start).count()
            << " ns" << endl;
        cout << chrono::duration_cast<chrono::microseconds>(end - start).count()
            << " ns" << endl;
        
        cout<< "K vecinos (digita K): ";
        cin >> k;
    }
    
    
    

    return 0;
}

