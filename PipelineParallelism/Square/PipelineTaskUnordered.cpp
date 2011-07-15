#include <iostream>
#include <cstdlib>
#include <ctime>

#include <omp.h>

#define TOKENS 50
#define WAITTIME_SEC 1

using namespace std;

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
  #pragma omp critical
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
  #pragma omp task untied
  {
    state = new Print(squared);
    state->next();
    delete state;
  }
  return NULL;
}
};

class Generate : public State {
private:
static int internalCounter;
public:
Generate() {
  #pragma omp critical
  {
    if(internalCounter < TOKENS) {
      internalCounter++;
      token = internalCounter;
      hasContents = true;
    }
  }
}

State* next() {
  State* state;

  int originalValue = token;
  #pragma omp task untied
  {
    // Original statement was
    // State* state = new Transform(token);
    // which caused an implicit reference to the this->token variable and the this pointer
    // The use of the this pointer seems to have caused some segmentation problems.
    state = new Transform(originalValue);
    state->next();
    delete state;
  }
  return NULL;
}
};

// Allocate a counter to count how many counters have been defined
// Modify/remove this as necesasry to facilitate parallelizing if it causes problems
int Generate::internalCounter = 0;


int main(int argc, char* argv[]) {
  srand( time(NULL) );

  cout << "Beginning simplest pipeline with fixed number of tokens" << endl;
  // We have fixed the number of tokens here but don't make the invalid assumption
  // that we will always know this in advance

  #pragma omp parallel
  {
    #pragma omp single
    {
      while(1) {
        State* state = new Generate();
        if( state->isValid() ) {
        #pragma omp task untied
          {
            state->next();
            delete state;
          }
        } else {
          delete state;
          break;
        }
      }
    }
  }

  return 0;
}
