//
// Created by shahnoor on 2/15/21.
//

#ifndef SQLATTICEPERCOLATION_ELEMENT_H
#define SQLATTICEPERCOLATION_ELEMENT_H


class Element_v13{
    int id = -1;
    int g_id = -1;

public:
    void reset(){
        g_id = -1;
//        id = -1;
    }
    void set_id(int id){
        this->id = id;
    }

    void set_gid(int gid){
        this->g_id = gid;
    }

    virtual int get_id(){
        return id;
    }

    virtual int get_gid(){
        return g_id;
    }


};


#endif //SQLATTICEPERCOLATION_ELEMENT_H
