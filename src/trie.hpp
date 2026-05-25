#include <unordered_map>
#include <vector>

using namespace std;

struct TrieNode {
    unordered_map<int, TrieNode*> next;
};

class Trie {
private:
    TrieNode* root = new TrieNode();

public:
    void insert(const vector<int>& arr) {
        TrieNode* node = root;

        for (int x : arr) {
            if (!node->next.count(x)) {
                node->next[x] = new TrieNode();
            }

            node = node->next[x];
        }
    }

    bool prefixExists(const vector<int>& prefix) {
        TrieNode* node = root;

        for (int x : prefix) {
            if (!node->next.count(x)) {
                return false;
            }

            node = node->next[x];
        }

        return true;
    }
};