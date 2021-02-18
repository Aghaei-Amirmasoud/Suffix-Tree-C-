// A C++ program to implement Ukkonen's Suffix Tree Construction
// Here we build generalized suffix tree for given string S
// and it's reverse R, then we find
// longest palindromic substring of given string S
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <unordered_set>
#include <fstream>

#define MAX_CHAR 256
using namespace std;

string input_1;
string input_2;
int len[2];

struct SuffixTreeNode {
    struct SuffixTreeNode *children[MAX_CHAR];

    //pointer to other node via suffix link
    struct SuffixTreeNode *suffixLink;

    /*(start, end) interval specifies the edge, by which the
    node is connected to its parent node. Each edge will
    connect two nodes, one parent and one child, and
    (start, end) interval of a given edge will be stored
    in the child node. Lets say there are two nods A and B
    connected by an edge with indices (5, 8) then this
    indices (5, 8) will be stored in node B. */
    int start;
    int *end;

    /*for leaf nodes, it stores the index of suffix for
    the path from root to leaf*/
    int suffixIndex;

    //To store indices of children suffixes in given string
    unordered_set<int> *forwardIndices;

    //To store indices of children suffixes in reversed string
    unordered_set<int> *reverseIndices;
};


typedef struct SuffixTreeNode Node;

char text[5368709]; //Input string
Node *root = nullptr; //Pointer to root node

/*lastNewNode will point to newly created internal node,
waiting for it's suffix link to be set, which might get
a new suffix link (other than root) in next extension of
same phase. lastNewNode will be set to NULL when last
newly created internal node (if there is any) got it's
suffix link reset to new internal node created in next
extension of same phase. */
Node *lastNewNode = nullptr;
Node *activeNode = nullptr;

/*activeEdge is represeted as input string character
index (not the character itself)*/
int activeEdge = -1;
int activeLength = 0;

// remainingSuffixCount tells how many suffixes yet to
// be added in tree
int remainingSuffixCount = 0;
int leafEnd = -1;
int *rootEnd = nullptr;
int *splitEnd = nullptr;
int size = -1;
int size1 = 0; //Size of 1st string
int reverseIndex; //Index of a suffix in reversed string
unordered_set<int>::iterator forwardIndex;

Node *newNode(int start, int *end)
{
    Node *node =(Node*) malloc(sizeof(Node));
    int i;
    for (i = 0; i < MAX_CHAR; i++)
        node->children[i] = nullptr;

    /*For root node, suffixLink will be set to NULL
    For internal nodes, suffixLink will be set to root
    by default in current extension and may change in
    next extension*/
    node->suffixLink = root;
    node->start = start;
    node->end = end;

    /*suffixIndex will be set to -1 by default and
    actual suffix index will be set later for leaves
    at the end of all phases*/
    node->suffixIndex = -1;
    node->forwardIndices = new unordered_set<int>;
    node->reverseIndices = new unordered_set<int>;
    return node;
}

int edgeLength(Node *n) {
    return *(n->end) - (n->start) + 1;
}

int walkDown(Node *currNode)
{
    /*activePoint change for walk down (APCFWD) using
    Skip/Count Trick (Trick 1). If activeLength is greater
    than current edge length, set next internal node as
    activeNode and adjust activeEdge and activeLength
    accordingly to represent same activePoint*/
    if (activeLength >= edgeLength(currNode))
    {
        activeEdge += edgeLength(currNode);
        activeLength -= edgeLength(currNode);
        activeNode = currNode;
        return 1;
    }
    return 0;
}

void extendSuffixTree(int pos)
{
    /*Extension Rule 1, this takes care of extending all
    leaves created so far in tree*/
    leafEnd = pos;

    /*Increment remainingSuffixCount indicating that a
    new suffix added to the list of suffixes yet to be
    added in tree*/
    remainingSuffixCount++;

    /*set lastNewNode to NULL while starting a new phase,
    indicating there is no internal node waiting for
    it's suffix link reset in current phase*/
    lastNewNode = nullptr;

    //Add all suffixes (yet to be added) one by one in tree
    while(remainingSuffixCount > 0) {

        if (activeLength == 0)
            activeEdge = pos; //APCFALZ

        // There is no outgoing edge starting with
        // activeEdge from activeNode
        if (activeNode->children[int(text[activeEdge])] == nullptr)
        {
            //Extension Rule 2 (A new leaf edge gets created)
            activeNode->children[int(text[activeEdge])] =
            newNode(pos, &leafEnd);

            /*A new leaf edge is created in above line starting
            from an existng node (the current activeNode), and
            if there is any internal node waiting for it's suffix
            link get reset, point the suffix link from that last
            internal node to current activeNode. Then set lastNewNode
            to NULL indicating no more node waiting for suffix link
            reset.*/
            if (lastNewNode != nullptr)
            {
                lastNewNode->suffixLink = activeNode;
                lastNewNode = nullptr;
            }
        }
        // There is an outgoing edge starting with activeEdge
        // from activeNode
        else
        {
            // Get the next node at the end of edge starting
            // with activeEdge
            Node *next = activeNode->children[int(text[activeEdge])] ;
            if (walkDown(next))//Do walkdown
            {
                //Start from next node (the new activeNode)
                continue;
            }
            /*Extension Rule 3 (current character being processed
            is already on the edge)*/
            if (text[next->start + activeLength] == text[pos])
            {
                if (lastNewNode != nullptr && activeNode != root)
                {
                    lastNewNode->suffixLink = activeNode;
                    lastNewNode = nullptr;
                }
                //APCFER3
                activeLength++;
                /*STOP all further processing in this phase
                and move on to next phase*/
                break;
            }

            /*We will be here when activePoint is in middle of
            the edge being traversed and current character
            being processed is not on the edge (we fall off
            the tree). In this case, we add a new internal node
            and a new leaf edge going out of that new node. This
            is Extension Rule 2, where a new leaf edge and a new
            internal node get created*/
            splitEnd = (int*) malloc(sizeof(int));
            *splitEnd = next->start + activeLength - 1;

            //New internal node
            Node *split = newNode(next->start, splitEnd);
            activeNode->children[int(text[activeEdge])] = split;

            //New leaf coming out of new internal node
            split->children[int(text[pos])] = newNode(pos, &leafEnd);
            next->start += activeLength;
            split->children[int(text[next->start])] = next;

            /*We got a new internal node here. If there is any
            internal node created in last extensions of same
            phase which is still waiting for it's suffix link
            reset, do it now.*/
            if (lastNewNode != nullptr)
            {
                /*suffixLink of lastNewNode points to current newly
                created internal node*/
                lastNewNode->suffixLink = split;
            }

            /*Make the current newly created internal node waiting
            for it's suffix link reset (which is pointing to root
            at present). If we come across any other internal node
            (existing or newly created) in next extension of same
            phase, when a new leaf edge gets added (i.e. when
            Extension Rule 2 applies is any of the next extension
            of same phase) at that point, suffixLink of this node
            will point to that internal node.*/
            lastNewNode = split;
        }

        /* One suffix got added in tree, decrement the count of
        suffixes yet to be added.*/
        remainingSuffixCount--;
        if (activeNode == root && activeLength > 0) //APCFER2C1
        {
            activeLength--;
            activeEdge = pos - remainingSuffixCount + 1;
        }
        else if (activeNode != root) //APCFER2C2
        {
            activeNode = activeNode->suffixLink;
        }
    }
}

void print(int i, int j)
{
    int k;
    for (k=i; k<=j && text[k] != '#'; k++) {
        printf("%c", text[k]);
    }
    if(k<=j) {
        printf("#");
    }
}

//Print the suffix tree as well along with setting suffix index
//So tree will be printed in DFS manner
//Each edge along with it's suffix index will be printed
void setSuffixIndexByDFS(Node *n, int labelHeight, const string& attach)
{
    if (n == nullptr) return;

    if (n->start != -1) //A non-root node
    {
        //Print the label on edge from parent to current node
        cout << attach.c_str();
        print(n->start, *(n->end));
    }
    int leaf = 1;
    int i;
    for (i = 0; i < MAX_CHAR; i++)
    {
        if (n->children[i] != nullptr)
        {

            if (leaf == 1 && n->start != -1)
                printf(" [%d]\n", n->suffixIndex);

            //Current node is not a leaf as it has outgoing
            //edges from it.
            leaf = 0;
            setSuffixIndexByDFS(n->children[i], labelHeight +
                                                edgeLength(n->children[i]), string(attach+"---"));
            if(n != root)
            {
                //Add chldren's suffix indices in parent
                n->forwardIndices->insert(
                        n->children[i]->forwardIndices->begin(),
                        n->children[i]->forwardIndices->end());
                n->reverseIndices->insert(
                        n->children[i]->reverseIndices->begin(),
                        n->children[i]->reverseIndices->end());
            }
        }
    }
    if (leaf == 1)
    {
        for(i= n->start; i<= *(n->end); i++)
        {
            if(text[i] == '#')
            {
                n->end = (int*) malloc(sizeof(int));
                *(n->end) = i;
            }
        }
        n->suffixIndex = size - labelHeight;
        if(n->suffixIndex < size1)  //Suffix of Given String
            n->forwardIndices->insert(n->suffixIndex);
        else //Suffix of Reversed String
            n->reverseIndices->insert(n->suffixIndex - size1);
            printf(" [%d]\n", n->suffixIndex);
    }
}

void freeSuffixTreeByPostOrder(Node *n)
{
    if (n == nullptr)
        return;
    int i;
    for (i = 0; i < MAX_CHAR; i++)
    {
        if (n->children[i] != nullptr)
        {
            freeSuffixTreeByPostOrder(n->children[i]);
        }
    }
    if (n->suffixIndex == -1) {
        free(n->end);
    }
    free(n);
}

/*Build the suffix tree and print the edge labels along with
suffixIndex. suffixIndex for leaf edges will be >= 0 and
for non-leaf edges will be -1*/
void buildSuffixTree()
{
    size = strlen(text);
    int i;
    rootEnd = (int*) malloc(sizeof(int));
    *rootEnd = - 1;

    /*Root is a special node with start and end indices as -1,
    as it has no parent from where an edge comes to root*/
    root = newNode(-1, rootEnd);

    activeNode = root; //First activeNode will be root
    for (i=0; i<size; i++)
        extendSuffixTree(i);
    int labelHeight = 0;
    setSuffixIndexByDFS(root, labelHeight, "");
}

void doTraversal_LPS(Node *n, int labelHeight, int* maxHeight,
                 int* substringStartIndex)
{
    if(n == nullptr)
    {
        return;
    }
    int i;
    if(n->suffixIndex < 0) //If it is internal node
    {
        for (i = 0; i < MAX_CHAR; i++)
        {
            if(n->children[i] != nullptr)
            {
                doTraversal_LPS(n->children[i], labelHeight +
                                            edgeLength(n->children[i]),
                            maxHeight, substringStartIndex);

                if(*maxHeight < labelHeight
                   && !n->forwardIndices->empty() &&
                   !n->reverseIndices->empty())
                {
                    for (forwardIndex=n->forwardIndices->begin();
                         forwardIndex!=n->forwardIndices->end();
                         ++forwardIndex)
                    {
                        reverseIndex = (size1 - 2) -
                                       (*forwardIndex + labelHeight - 1);
                        //If reverse suffix comes from
                        //SAME position in given string
                        //Keep track of deepest node
                        if(n->reverseIndices->find(reverseIndex) !=
                           n->reverseIndices->end())
                        {
                            *maxHeight = labelHeight;
                            *substringStartIndex = *(n->end) -
                                                   labelHeight + 1;
                            break;
                        }
                    }
                }
            }
        }
    }
}

void getLongestPalindromicSubstring()
{
    int maxHeight = 0;
    int substringStartIndex = 0;
    doTraversal_LPS(root, 0, &maxHeight, &substringStartIndex);

    int k;
    for (k=0; k<maxHeight; k++)
        printf("%c", text[k + substringStartIndex]);
    if(k == 0)
        printf("No palindromic substring");
    else
        printf(", of length: %d",maxHeight);
    printf("\n");
}

int doTraversal_LCS(Node *n, int labelHeight, int* maxHeight,
                int* substringStartIndex)
{
    if(n == nullptr)
    {
        return -1;
    }
    int i;
    int ret;
    if(n->suffixIndex < 0) //If it is internal node
    {
        for (i = 0; i < MAX_CHAR; i++)
        {
            if(n->children[i] != nullptr)
            {
                ret = doTraversal_LCS(n->children[i], labelHeight +
                                                  edgeLength(n->children[i]),
                                  maxHeight, substringStartIndex);

                if(n->suffixIndex == -1)
                    n->suffixIndex = ret;
                else if((n->suffixIndex == -2 && ret == -3) ||
                        (n->suffixIndex == -3 && ret == -2) ||
                        n->suffixIndex == -4)
                {
                    n->suffixIndex = -4;//Mark node as XY
                    //Keep track of deepest node
                    if(*maxHeight < labelHeight)
                    {
                        *maxHeight = labelHeight;
                        *substringStartIndex = *(n->end) -
                                               labelHeight + 1;
                    }
                }
            }
        }
    }
    else if(n->suffixIndex > -1 && n->suffixIndex < size1)//suffix of X
        return -2;//Mark node as X
    else if(n->suffixIndex >= size1)//suffix of Y
        return -3;//Mark node as Y
    return n->suffixIndex;
}

void getLongestCommonSubstring()
{
    int maxHeight = 0;
    int substringStartIndex = 0;
    doTraversal_LCS(root, 0, &maxHeight, &substringStartIndex);

    int k;
    for (k=0; k<maxHeight; k++)
        printf("%c", text[k + substringStartIndex]);
    if(k == 0)
        printf("No common substring");
    else
        printf(", of length: %d",maxHeight);
    printf("\n");
}

void doTraversal_LRS(Node *n, int labelHeight, int* maxHeight,
                 int* substringStartIndex)
{
    if(n == nullptr)
    {
        return;
    }
    int i;
    if(n->suffixIndex == -1) //If it is internal node
    {
        for (i = 0; i < MAX_CHAR; i++)
        {
            if(n->children[i] != nullptr)
            {
                doTraversal_LRS(n->children[i], labelHeight +
                                            edgeLength(n->children[i]), maxHeight,
                            substringStartIndex);
            }
        }
    }
    else if(n->suffixIndex > -1 &&
            (*maxHeight < labelHeight - edgeLength(n)))
    {
        *maxHeight = labelHeight - edgeLength(n);
        *substringStartIndex = n->suffixIndex;
    }
}

void getLongestRepeatedSubstring()
{
    int maxHeight = 0;
    int substringStartIndex = 0;
    doTraversal_LRS(root, 0, &maxHeight, &substringStartIndex);
    //printf("maxHeight %d, substringStartIndex %d\n", maxHeight,substringStartIndex);
    printf("Longest Repeated Substring in %s is: ", text);
    int k;
    for (k=0; k<maxHeight; k++)
        printf("%c", text[k + substringStartIndex]);
    if(k == 0)
        printf("No repeated substring");
    printf("\n");
}


int traverseEdge(const char *str, int idx, int start, int end)
{
    int k;
    //Traverse the edge with character by character matching
    for(k=start; k<=end && str[idx] != '\0'; k++, idx++)
    {
        if(text[k] != str[idx])
            return -1;  // mo match
    }
    if(str[idx] == '\0')
        return 1;  // match
    return 0;  // more characters yet to match
}

int doTraversalToCountLeaf(Node *n)
{
    if(n == nullptr)
        return 0;
    if(n->suffixIndex > -1)
    {
        printf("\nFound at position: %d", n->suffixIndex);
        return 1;
    }
    int count = 0;
    int i;
    for (i = 0; i < MAX_CHAR; i++)
    {
        if(n->children[i] != nullptr)
        {
            count += doTraversalToCountLeaf(n->children[i]);
        }
    }
    return count;
}

int countLeaf(Node *n)
{
    if(n == nullptr)
        return 0;
    return doTraversalToCountLeaf(n);
}

int doTraversal_PM(Node *n, char* str, int idx)
{
    if(n == nullptr)
    {
        return -1; // no match
    }
    int res;
    //If node n is not root node, then traverse edge
    //from node n's parent to node n.
    if(n->start != -1)
    {
        res = traverseEdge(str, idx, n->start, *(n->end));
        if(res == -1)  //no match
            return -1;
        if(res == 1) //match
        {
            if(n->suffixIndex > -1)
                printf("\nsubstring count: 1 and position: %d",
                       n->suffixIndex);
            else
                printf("\nsubstring count: %d", countLeaf(n));
            return 1;
        }
    }
    //Get the character index to search
    idx = idx + edgeLength(n);
    //If there is an edge from node n going out
    //with current character str[idx], traverse that edge
    if(n->children[str[idx]] != nullptr)
        return doTraversal_PM(n->children[str[idx]], str, idx);
    else
        return -1;  // no match
}

void checkForSubString(char* str)
{
    int res = doTraversal_PM(root, str, 0);
    if(res == 1)
        printf("\nPattern <%s> is a Substring\n", str);
    else
        printf("\nPattern <%s> is NOT a Substring\n", str);
}

void read_file_PM_LCS() {
    input_1 = "";
    fstream newfile;
    newfile.open("../IO/input_1.txt",ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()){   //checking whether the file is open
        string tp;
        int i = 0;
        while(getline(newfile, tp)){ //read data from file object and put it into string.
            if (tp != ">") {
                input_1 += tp + "#";
            }
            len[i] = tp.size();
            i++;
        }
        newfile.close(); //close the file object.
    }
    input_1.pop_back();
    input_1 += "$";
}

void read_file_LRS_LPS() {
    input_2 = "";
    fstream newfile;
    newfile.open("../IO/input_2.txt",ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()){   //checking whether the file is open
        string tp;
        while(getline(newfile, tp)){ //read data from file object and put it into string.
            if (tp != ">") {
                input_2 += tp;
            }
        }
        newfile.close(); //close the file object.
    }
    input_2 += "$";
}

// driver program to test above functions
int main(int argc, char *argv[])
{
    string tmp = "y";
    int choose;
    cin >> choose;
    if (choose == 1) {
        read_file_LRS_LPS();
        size1 = input_2.size();
        input_2.pop_back();
        cout << size1 << endl;
        string rev = string(input_2.rbegin(), input_2.rend());
        input_2 += "#" + rev + "$";
        cout << input_2 << endl;
        strcpy(text, input_2.c_str());
        buildSuffixTree();
        cout << endl;
        printf("Longest Palindromic Substring is: ");
        getLongestPalindromicSubstring();
        freeSuffixTreeByPostOrder(root);
    } else if (choose == 2 ) {
        freeSuffixTreeByPostOrder(root);
        read_file_LRS_LPS();
        strcpy(text, input_2.c_str());
        cout << endl;
        buildSuffixTree();
        cout << endl;
        getLongestRepeatedSubstring();
        freeSuffixTreeByPostOrder(root);
    } else if (choose == 3) {
        read_file_PM_LCS();
        strcpy(text, input_1.c_str());
        cout << endl;
        buildSuffixTree();
        cout << endl;
        size1 = len[1];
        cout << len[1] << endl;
        cout << input_1 << endl;
        printf("Longest Common Substring is: ");
        getLongestCommonSubstring();
        freeSuffixTreeByPostOrder(root);
    } else if (choose == 4) {
        freeSuffixTreeByPostOrder(root);
        read_file_PM_LCS();
        strcpy(text, input_1.c_str());
        cout << endl;
        buildSuffixTree();
        cout << endl;
        char pattern[] = "AT";
        checkForSubString(pattern);
        freeSuffixTreeByPostOrder(root);
    }

    return 0;
}




