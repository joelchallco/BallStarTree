#ifndef __NODE_BALL_STAR_H
#define __NODE_BALL_STAR_H
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

struct rid {
public:
	int page;
	int slot;

	rid(int _page = -1,int _slot = -1) {
		page = _page;
		slot = _slot;
		}
	
	void set_rid(int _page,int _slot) {
		page = _page;
		slot = _slot;
		}

};


class NodeBallStar{
public:
    vector<float> data;
	string name;
    float center;
    float radio;

    NodeBallStar* left;
    NodeBallStar* right;

    NodeBallStar();
    void getcenter();
    void getradio();

};




#endif