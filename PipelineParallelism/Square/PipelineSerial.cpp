#include <iostream>
#include <cstdlib>
#include <ctime>

#include <omp.h>

#define TOKENS 50

using namespace std;

int WAITTIME_SEC = 1;

void wait(int seconds) {
  double endwait = omp_get_wtime() + seconds;

  while(omp_get_wtime() < endwait) { }
}

class State {
protected:
int token;
bool hasContents;
public:
virtual State* next() = 0;
bool isValid() {
  return hasContents;
}

State () {
	token = -1;
  hasContents = false;
}
virtual ~State() { }
};

class Print : public State {
public:
Print(int _token) {
  token = _token;
}

State* next() {
  cout << "(" << token << ")" << endl;
  return NULL;
}
};

class Transform : public State {
public:
Transform(int _token) {
  token = _token;
}
State* next() {
  int squared = token * token;

  wait(WAITTIME_SEC);

  State* state;
  state = new Print(squared);
  state->next();
  delete state;
  return NULL;
}
};

class Generate : public State {
private:
static int internalCounter;
public:
Generate() {
  if(internalCounter < TOKENS) {
    internalCounter++;
    token = internalCounter;
    hasContents = true;
  }
}

State* next() {
  State* state = new Transform(token);

  state->next();
  delete state;
  return NULL;
}
};

// Allocate a counter to count how many counters have been defined
// Modify/remove this as necesasry to facilitate parallelizing if it causes problems
int Generate::internalCounter = 0;


int main(int argc, char* argv[]) {
  srand( time(NULL) );

  if(argc >= 2)
    WAITTIME_SEC = atoi(argv[1]);

  cout << "Beginning simplest pipeline with fixed number of tokens" << endl;
  cout << "Transform stage takes: " << WAITTIME_SEC << endl;
  // We have fixed the number of tokens here but don't make the invalid assumption
  // that we will always know this in advance

  while(1) {
    State* state = new Generate();
    if( state->isValid() ) {
      state->next();
      delete state;
    } else {
      delete state;
      break;
    }
  }

  return 0;
}
