#include "MyIGtrie.h"

MyIGtrie::MyIGtrie()
{

}

MyIGtrie::~MyIGtrie()
{
  destroy();
}

void MyIGtrie::init(int K)
{
  maxLabels = 10 * K * K;
  numLabels = 1;
  labelPaths = (int**) malloc(sizeof(int*) * maxLabels);
  labelLeaf = (int*) malloc(sizeof(int) * maxLabels);
  labelCount = (int*) malloc(sizeof(int) * maxLabels);

  labelPaths[0] = new int[LB_WORD_SIZE];
  labelLeaf[0] = 1;
  labelCount[0] = 0;
  memset(labelPaths[0], -1, sizeof(int) * LB_WORD_SIZE);
}

void MyIGtrie::destroy()
{
  if (!numLabels)
    return;

  int i;
  for (i = 0; i < numLabels; i++)
    delete[] labelPaths[i];

  free(labelPaths);
  free(labelLeaf);
  free(labelCount);
  numLabels = 0;
  maxLabels = 0;
}

void MyIGtrie::expand()
{
  maxLabels *= 2;
  labelPaths = (int**) realloc(labelPaths, sizeof(int*) * maxLabels);
  labelLeaf = (int*) realloc(labelLeaf, sizeof(int) * maxLabels);
  labelCount = (int*) realloc(labelCount, sizeof(int) * maxLabels);
}

void MyIGtrie::incrementLabel(int labelNode, int value)
{
  labelCount[labelNode] += value;
}

int MyIGtrie::insertLabel(int labelNode, long long int label, int digits, bool createNew)
{
  if (!createNew) {
    return labelPaths[labelNode][label & (LB_WORD_SIZE - 1)];
  }

  if (labelPaths[labelNode][label & (LB_WORD_SIZE - 1)] == -1)
  {
    if (numLabels == maxLabels)
      expand();

    int newNode = numLabels++;
    labelPaths[newNode] = new int[LB_WORD_SIZE];
    memset(labelPaths[newNode], -1, sizeof(int) * LB_WORD_SIZE);
    labelPaths[labelNode][label & (LB_WORD_SIZE - 1)] = newNode;
    labelLeaf[newNode] = ((digits <= LB_WORD_LEN) ? digits : 0);
    labelCount[newNode] = 0;
  }

  int nextNode = labelPaths[labelNode][label & (LB_WORD_SIZE - 1)];

  if (!labelLeaf[nextNode])
    return insertLabel(nextNode, label >> LB_WORD_LEN, digits - LB_WORD_LEN);
  return nextNode;
}

vector<pair<long long int, int> > MyIGtrie::enumerate(int K)
{
  enumeration.clear();
  enumerateFrom(0, 0, 0, 0, K - 1);

  return enumeration;
}

void MyIGtrie::enumerateFrom(int currentNode, long long int label, long long int parLabel, int parSize, int remaining)
{
  if (remaining == 0)
  {
    enumeration.push_back(make_pair(label, labelCount[currentNode]));
    return;
  }

  int i;
  for (i = 0; i < LB_WORD_SIZE; i++)
    if (labelPaths[currentNode][i] != -1)
    {
      long long int tmpLabel = parLabel;
      int tmpSize = parSize;

      if (labelLeaf[labelPaths[currentNode][i]])
      {
        int digits = labelLeaf[labelPaths[currentNode][i]];
        tmpLabel |= (i << tmpSize);
        tmpSize += digits;
        enumerateFrom(labelPaths[currentNode][i], ((label << tmpSize) | tmpLabel), 0, 0, remaining - 1);
      }
      else
      {
        int digits = LB_WORD_LEN;
        tmpLabel |= (i << tmpSize);
        tmpSize += digits;
        enumerateFrom(labelPaths[currentNode][i], label, tmpLabel, tmpSize, remaining);
      }
    }
}
