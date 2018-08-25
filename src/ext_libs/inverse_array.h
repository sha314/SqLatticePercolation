//
// Created by shahnoor on 11/23/2017.
//

#ifndef TREES_INVERSE_ARRAY_H
#define TREES_INVERSE_ARRAY_H


#include <cstdlib>
#include <iostream>
#include <map>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

/**
 *
 * todo Implement Iterator
 */

/**
 * Node for double linked list for inverse array
 */
template <typename T>
struct NodeDIA {
    T key_;
    size_t order_;
    NodeDIA *prev_{nullptr}; // pointer to previous node
    NodeDIA *next_{nullptr}; // pointer to next node

    ~NodeDIA() = default;
//    size_t order() const {return order_;}
};



/*************************************************************************************
 * Basically a Doublylinked list
 * But does not have insert method
 * Instead values are always numbers starting from 0
 * If a node is deleted values are reordered from that point
 *
 * Designed to be used in InverseArray class
 *
 * organise inserted keys in the order of insersion,
 * first value is indexed as 0 and second
 * value as 1 and so on. when called using
 * value it returns the index of that value.
 * Inverse senario of an array.
 *
 * Insert operation for key only. value will be set automatically
 * in the order of insertion starting from 0.
 *
 * Erase operation will erase one key and the value associated with it.
 * All the value associated with the later key sequence will decrease by 1.
 *
 * e.g., using  key->value convention
 *
 * before :
 * 20->0
 * 45->1
 * 5->2
 * 10->3
 * 50->4
 *
 * erase key 5
 * after:
 * 20->0
 * 45->1
 * 10->2
 * 50->3
 *
 */
template <typename T>
class LinkedIndex {
private:
    size_t counter = {};    // counts number of elements present in the object
    NodeDIA<T> *head, *tail;
    NodeDIA<T> * iterator;

public:

    ~LinkedIndex()
    {
        clear();
    }

    LinkedIndex() {
        head = nullptr;
        tail = nullptr;
    }

//    NodeDIA<T>& operator=(const NodeDIA<T>& other) = delete;
//    NodeDIA<T>& operator=(NodeDIA<T>&& other) = delete;

    void clear();

    /**************
     * Iterator
     *************/
//    NodeDIA<T>* begin(){
//        iterator = head;
//        return iterator;
//    }
//
//    NodeDIA<T>* end(){
//        return tail;
//    }
//
//    NodeDIA<T>* operator++(){
//        if(iterator >= tail){
//            return new NodeDIA<T>;
//        }
//        iterator = iterator->next_;
//        return iterator;
//    }
//
//    bool operator!=(const NodeDIA<T>& itr) const {
//        return *this != itr;
//    }

    /**************************************
     * Create new node at the end of the list
     */

    NodeDIA<T>* create_node(size_t key);

    /***********************************
     * Display
     ********************************/
    void display();

    /**
     * See how the Nodes are arranged.
     *
     */
    void view_fdw();        // Fowrward view
    void view_bwd();        // Backward view

    std::string str() const;    // InverseArray as string


    /****************************************
     * Insertion
     **************************************/

    NodeDIA<T>* insert_start(size_t value);
    NodeDIA<T>* insert_last(size_t value);
    NodeDIA<T>* insert_position(size_t position, size_t value);

    /****************************************
     * Erasing and Deletion
     ***************************************/
    size_t erase(NodeDIA<T>* n);
    void delete_first();
    void delete_last();
    void delete_position(size_t position);
    void delete_key(size_t key);    // Find and delete given key


    /************************************************
     * Swapping
     ***********************************************/

    /**
     * Insert n2 in the position n1 and then erase n2
     * Condition: n1 and n2 must exists in the LinkedIndex object
     * @param n1
     * @param n2
     */
    void insert_in_place(NodeDIA<T>* n1, NodeDIA<T>* n2);


    /**
     * Insert a new Node which is n2 in position n1
     * n1 will be pushed downward
     * @param n1
     * @param n2
     */
    void insert_new_in_place(NodeDIA<T>* n1, NodeDIA<T>* n2);


    /**
     * Condition: both node must exists in the LinkedIndex object
     * @param n1
     * @param n2
     */
    void swap(NodeDIA<T>* n1, NodeDIA<T>* n2);


private:
    /**
     * Reorder from next node.
     * Given that current node is in order with its previous node
     * Must be called internally. Uncomment the commented code if neccessary.
     * @param n
     */
    void re_order(NodeDIA<T> *n);

};



/************************************************
 * LinkedIndex
 * method definitions
 ************************************************/

/**
 * todo // test more use internet
 */
template <typename T>
void LinkedIndex<T>::clear(){
    while(head != nullptr){
        NodeDIA<T>* temp = head;
        head = head->next_;
//            cout << "deleting " << temp->key << endl;
        delete temp;  // delete previous, i.e., which was current a step earlier
        --counter;
    }
//    tail = nullptr;
//    head = nullptr;

//    while (head->next_ != nullptr) {
//        NodeDIA<T>* temp = head;
//        head = head->next_;
//        delete temp;
//    }

    tail = nullptr;
    head = nullptr;
}

/************
 * Create new node at the end of the list
 */
/**
 *
 * @param key
 * @return address where current key is placed
 */
template <typename T>
NodeDIA<T>* LinkedIndex<T>::create_node(size_t key) {
    NodeDIA<T> *temp = new NodeDIA<T>;
    temp->key_ = key;
    ++counter;
    if (head == nullptr) {
        temp->order_ = 0;    // initial point
        head = temp;
        tail = temp;
//            delete temp;
        return head;
    } else {
        temp->order_ = tail->order_ + 1;
        tail->next_ = temp;
        temp->prev_ = tail;
        tail = temp;
        return tail;
    }
}

/**********
 * Display
 */
template <typename T>
void LinkedIndex<T>::display()
{
    NodeDIA<T> *temp= new NodeDIA<T>;
    std::cout << "Total element " << counter << endl;
    std::cout << "Forward : " << std::endl;
    temp=head;
    size_t tmp_count{};
    while(temp != nullptr)
    {
        std::cout<<temp->key_<<", " << temp->order_ << endl;
        temp=temp->next_;
        ++tmp_count;
        if(tmp_count == counter){
//                cout << "Reached" << endl;
            break;
        }
    }
    cout << "Total element " << tmp_count << endl;

    std::cout << std::endl << "Backward : " << std::endl;
    temp=tail;
    tmp_count = 0;
    while(temp != nullptr)
    {
        std::cout<<temp->key_<<", " << temp->order_ << endl;
        temp=temp->prev_;
        ++tmp_count;
        if(tmp_count == counter){
//                cout << "Reached" << endl;
            break;
        }
    }
    cout << "Total element " << tmp_count << endl;
}

/**
 * See how the Nodes are arranged.
 * Fowrward view
 */
template <typename T>
void LinkedIndex<T>::view_fdw()
{
    NodeDIA<T> *temp = new NodeDIA<T>;
    temp=head;
    std::cout << "Forward viewing "<< endl;
    std::cout << "Total element " << counter << endl;
    std::cout << "previous address <= current key, order => next address" << endl;
    size_t tmp_count{};
    while(temp != nullptr)
    {
        std::cout << std::setw(10) << temp->prev_ <<" <prev=("
                  << std::setw(6) << temp->key_ << ", " << std::setw(6) << temp->order_
                  << ")=next> " << std::setw(10) << temp->next_ << std::endl;
        temp=temp->next_;
        ++tmp_count;
        if(tmp_count == counter){
            cout << "Reached" << endl;
            break;
        }
    }
    delete temp;
}

/**
 * Backward view
 */
template <typename T>
void LinkedIndex<T>::view_bwd()
{
    NodeDIA<T> *temp = new NodeDIA<T>;
    temp=tail;
    std::cout << "Backward viewing "<< endl;
    std::cout << "Total element " << counter << endl;
    std::cout << "previous address <= current key, order => next address" << endl;
    size_t tmp_count{};
    while(temp != nullptr)
    {
        std::cout << std::setw(10) << temp->prev_ <<" <prev=("
                  << std::setw(6) << temp->key_ << ", " << std::setw(6) << temp->order_
                  << ")=next> " << std::setw(10) << temp->next_ << std::endl;
        temp=temp->prev_;
        ++tmp_count;
        if(tmp_count == counter){
            cout << "Reached" << endl;
            break;
        }
    }
    delete temp;
}

/**
 *
 * @return
 */
template <typename T>
std::string LinkedIndex<T>::str() const {
    stringstream ss;
    NodeDIA<T> *temp= new NodeDIA<T>;
    std::cout << "Total element " << counter << endl;
    temp=head;
    size_t tmp_count{};
    while(temp != nullptr)
    {
        ss << temp->key_<<", " << temp->order_ << endl;
        temp=temp->next_;
        ++tmp_count;
        if(tmp_count == counter){
//                cout << "Reached" << endl;
            break;
        }
    }
    return ss.str();
}

/**************
 * Insertion
 */

/**
 * insert in the begining of the list
 * @param value
 * @return
 */
template <typename T>
NodeDIA<T>* LinkedIndex<T>::insert_start(size_t value)
{
    NodeDIA<T> *temp=new NodeDIA<T>;

    temp->key_=value;

    temp->next_ = head;
    head->prev_ = temp;

    head=temp;
    ++counter;
    re_order(head);
    return head;
}

/**
 * insert at the last position
 * @param value
 * @return
 */
template <typename T>
NodeDIA<T>* LinkedIndex<T>::insert_last(size_t value){
    NodeDIA<T>* temp = new NodeDIA<T>;
    temp->key_ = value;
    temp->order_ = tail->order_ + 1;

    tail->next_ = temp;
    temp->prev_ = tail;
    tail = temp;
    tail->next_ = nullptr;
    ++counter;
    re_order(tail->prev_);
    return tail;
}

/**
 * insert at a particular position
 * @param pos starts from 0
 * @param value
 */
template <typename T>
NodeDIA<T>* LinkedIndex<T>::insert_position(size_t pos, size_t value)
{
    if(pos == 0){
        return insert_start(value);
    }
    if(pos == counter){
        return insert_last(value);
    }
    if(pos > counter){
        std::cerr << "Outside the range : line " << __LINE__ << std::endl;
        return new NodeDIA<T>;
    }
    NodeDIA<T> *previous   =new NodeDIA<T>;
    NodeDIA<T> *current    =new NodeDIA<T>;
    NodeDIA<T> *temp       =new NodeDIA<T>;
    current=head;
    for(size_t i=0;i<pos;i++)
    {
        previous=current;
        current=current->next_;
    }
    temp->key_=value;

    previous->next_=temp;
    temp->next_=current;

    temp->prev_ = previous;
    current->prev_ = temp;

    ++counter;
    re_order(previous);
    return temp;
}

/*****
 * Erasing using pointer
 */

/**
 * Erase given node
 * Cannot do this if we cannot access the previous node
 * @param n -> given node
 * @return current order of the erasing key
 */
template <typename T>
size_t LinkedIndex<T>::erase(NodeDIA<T>* n){
//        std::cout << "erasing " << n->key_ << ", " << n->order_  << " : line " << __LINE__ << std::endl;
    size_t a = n->order_;
    if(n == head){
        NodeDIA<T>* next = head;   // n contains head
        head = head->next_;
        head->prev_ = nullptr;  // setting first element as nullptr
        delete next;
        --counter;
        head->order_ = 0;
        re_order(head);
        return a;
    }
    if(n == tail){
        NodeDIA<T>* previous = tail;
        tail = tail->prev_;
        tail->next_ = nullptr;
        delete previous;
        --counter;
        return a;
    }
    NodeDIA<T>* previous = n->prev_;
    NodeDIA<T>* next = n->next_;
    previous->next_ = next;
    next->prev_ = previous;
    re_order(previous);
    delete n;
    --counter;
    return a;
}

/*************************
 * Deletion
 ***********************/
/**
 * Delete first element
 */
template <typename T>
void LinkedIndex<T>::delete_first()
{
    NodeDIA<T> *temp=new NodeDIA<T>;
    temp=head;
    head=head->next_;
    head->prev_ = nullptr;
    delete temp;
    --counter;
    head->order_ = 0;
    re_order(head);
}

/**
 * Delete last element
 */
template <typename T>
void LinkedIndex<T>::delete_last()
{
    NodeDIA<T> *current=new NodeDIA<T>;
    current=tail;
    tail = tail->prev_;
    tail->next_ = nullptr;
    delete current;
    --counter;
}


/**
 * Delete element in specified position and then reorder them
 * @param pos starts from 0
 */
template <typename T>
void LinkedIndex<T>::delete_position(size_t pos)
{
    NodeDIA<T> *current =new NodeDIA<T>;
    NodeDIA<T> *previous=new NodeDIA<T>;
    current=head;
    for(size_t i=0;i<pos;i++)
    {
        previous=current;
        current=current->next_;
    }
    previous->next_ = current->next_;
    current->prev_ = previous->prev_;
    --counter;
    delete current;
    re_order(previous);
}

/**
 * Find and delete given key
 * @param key
 */
template <typename T>
void LinkedIndex<T>::delete_key(size_t key)
{
    NodeDIA<T> *previous = nullptr;
    NodeDIA<T> *temp = nullptr;
    NodeDIA<T>* next = nullptr;
    previous=head;
    while(previous != nullptr)
    {
        if (previous->next_->key_ == key){
            // key found
            temp = previous->next_;
            next = temp->next_;
            previous->next_ = next;
            next->prev_ = previous;
            delete temp;
            re_order(previous);
            break;
        }
        previous=previous->next_;
    }
}

/**
 * Swapping
 */


/**
 * Insert n2 in the position n1 and then erase n2
 * Condition: n1 and n2 must exists in the LinkedIndex object
 * @param n1
 * @param n2
 */
template <typename T>
void LinkedIndex<T>::insert_in_place(NodeDIA<T>* n1, NodeDIA<T>* n2){
    n1->key_ = n2->key_;

    if(n1->order_ < n2->order_) {
        re_order(n1->prev_);
    }
    else{
        re_order(n2->prev_);
    }
    erase(n2);

}

/**
 * Insert a new Node which is n2 in position n1
 * n1 will be pushed downward
 * @param n1
 * @param n2
 */
template <typename T>
void LinkedIndex<T>::insert_new_in_place(NodeDIA<T>* n1, NodeDIA<T>* n2){
    // . . . n1->prev n1 . . .
    // to
    // . . . n1->prev n2 n1 . . .
    // place n2 in between n1->prev and n1
    NodeDIA<T>* a = n1->prev_;
    NodeDIA<T>* b = n1;    // copy of n1

    a->next_ = n2;
    b->prev_ = n2;
    n2->next_ = b;
    n2->prev_ = a;

    ++counter;
    re_order(a);
}

/**
 * Condition: both node must exists in the LinkedIndex object
 * @param n1
 * @param n2
 */
template <typename T>
void LinkedIndex<T>::swap(NodeDIA<T>* n1, NodeDIA<T>* n2){
    // just swap the keys
    size_t key1 = n1->key_;
    n1->key_ = n2->key_;
    n2->key_ = key1;
}

/**
 * Reorder from next node.
 * Given that current node is in order with its previous node
 * Must be called internally. Uncomment the commented code if neccessary.
 * @param n
 */
template <typename T>
void LinkedIndex<T>::re_order(NodeDIA<T> *n){
    size_t order = n->order_;
    while(n->next_ != nullptr){
        n->next_->order_ = ++order;
        n = n->next_;
//            if(order >= counter){
//                std::cout << "order >= number of elements" << endl;
//                break;
//            }
    }
}

/**
 *
 */
struct RetrieveKey
{
    template <typename T>
    typename T::first_type operator()(T keyValuePair) const
    {
        return keyValuePair.first;
    }
};


/******************************************************************************
 *
 * class: InverseArray
 * Used Algorithms: Binary Tree, Doubly Linked List
 * Binary Tree -> For fast searching
 * Doubly Linked List -> To order keys in the sequence of insertion
 *
 *
 * organise them in the order of insertion,
 * first value is indexed as 0 and second
 * value as 1 and so on. when called using
 * value it returns the index of that value.
 * Inverse senario of an array.
 *
 * Insert operation for key only. value will be set automatically
 * in the order of insertion starting from 0.
 *
 * Erase operation will erase one key and the value associated with it.
 * All the value associated with the later key sequence will decrease by 1.
 *
 * e.g., using  key->value convention
 *
 * before :
 * 20->0
 * 45->1
 * 5->2
 * 10->3
 * 50->4
 *
 * erase key 5
 * after:
 * 20->0
 * 45->1
 * 10->2
 * 50->3
 *
 *
 * todo // add iterator
 */
template <typename T>
class InverseArray {
    std::map<size_t, NodeDIA<T>*> _tracker;
    LinkedIndex<T> _ll;

public:
    ~InverseArray() = default;
    InverseArray() = default;


    /***********************************************
     * Iterator
     **********************************************/
//    iterator begin(){
//        return _tracker.begin();
//    }

    /**
     *
     */
    void clear()
    {
        _tracker.clear();
        _ll.clear();
    }

    /**
     *
     * @return
     */
    bool empty()
    {
        return _tracker.empty();
    }

    /**
     *
     * @return
     */
    size_t count(T k)
    {
        return _tracker.count(k);
    }

    /**
     *
     * @return -> array of keys
     */
    std::vector<T> keys(){
        vector<int> keys;

// Retrieve all keys
        std::transform(_tracker.begin(), _tracker.end(), back_inserter(keys), RetrieveKey());

// Dump all keys
//        std::copy(keys.begin(), keys.end(), std::ostream_iterator<int>(cout, "\n"));

        return keys;
    }

    /**
     *
     * @param key
     */
    void insert(T key);
//    {
//        if(_tracker.count(key) > 0){
//            // key already exists
//            std::cout << "Key already exits : line " << __LINE__ << std::endl;
//            return;
//        }
//        NodeDIA<T>* n = _ll.create_node(key);
//        _tracker[key] = n;
//    }

    /**
     * Insert key and order
     * and place key, order pair in correct place then
     * call re_order function but this should not change order that is given as argument
     * @param key
     * @param order
     */
    void insert_key_order(T key, size_t order){

    }

    /**
     *
     * @param key
     * @return
     */
    size_t find(T key){
        return _tracker[key]->order_;
    }

    /**
    *
    * @param key
    * @return
    */
    size_t operator[](T key) {
        // todo if key does not exists
        if( _tracker[key] == nullptr){
            std::cout << "if( _tracker[key] == nullptr) : line " << __LINE__ << std::endl;
            std::cout << "calling exit(1)" << endl;
            exit(1);
        }
        return _tracker[key]->order_;
    }

    /**
     *
     * @param key
     */
    size_t erase(T key){
        size_t a = _ll.erase(_tracker[key]);
        _tracker.erase(key);
        return a;
    }

    /**
     *
     * @param key1
     * @param key2
     */
    void swap(T key1, T key2){
        _ll.swap(_tracker[key1], _tracker[key2]);
    }

    /**
     *
     * @return
     */
    std::string str() const {
        return _ll.str();
    }

    /**
     *
     * @param os
     * @param x
     * @return
     */
    friend std::ostream& operator<<(std::ostream& os, const InverseArray<T> & x){
        return os << x.str() << endl;
    }
};

/************************************************
 * InverseArray
 * Methods definitions
 *************************************************/

/**
 *
 * @param key
 */
template <typename T>
void InverseArray<T>::insert(T key){
    if(_tracker.count(key) > 0){
        // key already exists
        std::cerr << "InverseArray<T>::insert -> Key already exits. line " << __LINE__ << std::endl;
        return;
    }
    NodeDIA<T>* n = _ll.create_node(key);
    _tracker[key] = n;
}


#endif //TREES_INVERSE_ARRAY_H
