#include <iostream>
#include <cstdlib>
#include <ctime>

#include <utility>
#include <queue>

#include <omp.h>

#define TOKENS 50
#define WAITTIME_SEC 1

using namespace std;

typedef pair<int,int> WrappedToken;

void wait(int seconds) {
  double endwait = omp_get_wtime() + seconds;

  while(omp_get_wtime() < endwait) { }
}

class State {
protected:
WrappedToken token;
bool hasContents;
public:
virtual State* next() = 0;
bool isValid() {
  return hasContents;
}
State () {
  token = make_pair(0, 0);
  hasContents = false;
}

virtual ~State() { }
};

// By default pair<T,T> are sorted in the priority queue such that the highest
// value (pair.first) occurs at the top of the queue. We need to reverse the order
class CustomPairComparer {
public:
bool operator()(const WrappedToken& lhs, const WrappedToken& rhs) {
  return lhs.first > rhs.first;
}
};

class Print : public State {
private:
static priority_queue<WrappedToken, vector<WrappedToken>, CustomPairComparer> internalBuffer;
static int tokenIDToPrint;
public:
Print(WrappedToken _token) {
  token = _token;
}

State* next() {
  #pragma omp critical
  {
    //Put this into the queue
    internalBuffer.push(token);
  }

  WrappedToken pairToPrint;
  bool shouldPrint = false;

  while(1) {
    #pragma omp critical
    {
      if( !internalBuffer.empty() ) {
        pairToPrint = internalBuffer.top();
        if(pairToPrint.first == tokenIDToPrint) {
          internalBuffer.pop();
          shouldPrint = true;
        }
      }
    }
    if(shouldPrint) {
      #pragma omp critical
      {
        cout << "(" << pairToPrint.first << "," << pairToPrint.second << ")" << endl;
        tokenIDToPrint++;
      }
      return NULL;
    }
  }
}
};

// Include a static queue for the printStage
priority_queue<WrappedToken, vector<WrappedToken>, CustomPairComparer> Print::internalBuffer;
// Include a static counter for the print stage to keep track of tokens
int Print::tokenIDToPrint = 1;

class Transform : public State {
public:
Transform(WrappedToken _token) {
  token = _token;
}
State* next() {
  int originalID = token.first;
  int squared = token.second * token.second;

  wait(WAITTIME_SEC);

  #pragma omp task untied
  {
    State* state = new Print( WrappedToken(originalID, squared) );
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
  if(internalCounter < TOKENS) {
    internalCounter++;
    token = make_pair(internalCounter, internalCounter);
    hasContents = true;
  }
}

State* next() {
  int originalID = token.first;
  int originalValue = token.second;

  #pragma omp task untied
  {
    // Original statement was
    // State* state = new Transform(token);
    // which caused an implicit reference to the this->token variable and the this pointer
    // The use of the this pointer seems to have caused some segmentation problems.
    State* state = new Transform( WrappedToken(originalID, originalValue) );
    state->next();
    delete state;
  }
  return NULL;
}
};

// Allocate a counter to count how many counters have been defined
// Modify/remove this as necessary to facilitate parallelizing if it causes problems
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
