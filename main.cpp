#include <iostream>
#include <regex>
#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>
using namespace std;
struct Right_Matrix{
    double value;
    bool is_node;
    int index;
};
struct Matrices{
    vector<vector<double>> Left={};
    vector<Right_Matrix> Right={};
    vector<vector<double>> L_Left={};
    vector<vector<double>> U_Left={};
    vector<double> Answer={};
};
class Matrix_solve{
public:
    Matrices Primary_TRAN;
    Matrices TRAN;
    Matrices DC;
    bool LUsetter(Matrices &M){
        int size=M.Left.size();
        for(int i=0;i<size;i++){
            bool All_are_zero=1;
            for(int j=0;j<size;j++) {
                if(M.Left[i][j]!=0){
                    All_are_zero=0;
                    break;
                }
            }
            if(All_are_zero) return false;
            All_are_zero=1;
            for(int j=0;j<size;j++) {
                if(M.Left[j][i]!=0){
                    All_are_zero=0;
                    break;
                }
            }
            if(All_are_zero) return false;
        }
        M.L_Left.resize(size);
        M.U_Left.resize(size);
        for(int i=0;i<size;i++){
            M.L_Left[i].clear();
            M.U_Left[i].clear();
            M.L_Left[i].resize(size,0.0);
            M.U_Left[i].resize(size,0.0);
            M.L_Left[i][i]=1.0;
        }
        for(int layer=0;layer<size;layer++){
            bool swap_is_okay= false;
            for(int swap_index=layer;swap_index<size;swap_index++){
                double diagonal_value=M.Left[swap_index][layer];
                for(int i=0;i<layer;i++) {
                    diagonal_value-=M.L_Left[swap_index][i]*M.U_Left[i][layer];
                }
                if(fabs(diagonal_value)>1e-12){
                    swap_is_okay= true;
                    swap(M.Left[layer],M.Left[swap_index]);
                    swap(M.Right[layer],M.Right[swap_index]);
                    swap(M.L_Left[layer],M.L_Left[swap_index]);
                    M.L_Left[layer][swap_index]=0.0;
                    M.L_Left[swap_index][layer]=0.0;
                    M.L_Left[swap_index][swap_index]=1.0;
                    M.L_Left[layer][layer]=1.0;
                    for(int half_layer=layer;half_layer<size;half_layer++){
                        double value=M.Left[layer][half_layer];
                        for(int i=0;i<layer;i++) value-=M.L_Left[layer][i]*M.U_Left[i][half_layer];
                        M.U_Left[layer][half_layer]=value;
                        if(half_layer!=layer){
                            value=M.Left[half_layer][layer];
                            for(int i=0;i<layer;i++) value-=M.L_Left[half_layer][i]*M.U_Left[i][layer];
                            M.L_Left[half_layer][layer]=value/M.U_Left[layer][layer];
                        }
                    }
                    break;
                }
            }
            if(!swap_is_okay) return false;
        }
        return true;
    }
    void Solve(Matrices &M){
        int size = M.Right.size();
        vector<double> external_answer(size);
        M.Answer.resize(size);
        for(int i=0;i<size;i++){
            external_answer[i]=M.Right[i].value;
            for(int j=0;j<i;j++){
                external_answer[i]-=M.L_Left[i][j]*external_answer[j];
            }
        }
        for(int i=size-1;i>=0;i--){
            M.Answer[i]=external_answer[i];
            for(int j=i+1;j<size;j++){
                M.Answer[i]-=M.U_Left[i][j]*M.Answer[j];
            }
            M.Answer[i]/=M.U_Left[i][i];
        }
    }
};
class Node{
public:
    string name;
    int index;
    bool is_ground=0;
    void Add_Equation(Matrix_solve &m){
        vector<double> add={};
        for(int i=0;i<m.Primary_TRAN.Left.size();i++) m.Primary_TRAN.Left[i].insert(m.Primary_TRAN.Left[i].begin()+index,0.0);
        add.resize(m.Primary_TRAN.Left.size()+1,0.0);
        m.Primary_TRAN.Left.push_back(add);
        m.Primary_TRAN.Right.push_back({0.0,1,index});
        for(int i=0;i<m.DC.Left.size();i++) m.DC.Left[i].insert(m.DC.Left[i].begin()+index,0.0);
        add.resize(m.DC.Left.size()+1,0.0);
        m.DC.Left.push_back(add);
        m.DC.Right.push_back({0.0,1,index});
    }
    void Add_GND(Matrix_solve &m){
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==index){
                for(int j=0;j<m.Primary_TRAN.Left[i].size();j++){
                    m.Primary_TRAN.Left[i][j]=0.0;
                }
                m.Primary_TRAN.Left[i][index]=1.0;
                m.Primary_TRAN.Right[i].value=0.0;
            }
        }
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==index){
                for(int j=0;j<m.TRAN.Left[i].size();j++){
                    m.TRAN.Left[i][j]=0.0;
                }
                m.TRAN.Left[i][index]=1.0;
                m.TRAN.Right[i].value=0.0;
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==index){
                for(int j=0;j<m.DC.Left[i].size();j++){
                    m.DC.Left[i][j]=0.0;
                }
                m.DC.Left[i][index]=1.0;
                m.DC.Right[i].value=0.0;
            }
        }
    }
    void Delete_Equation(Matrix_solve &m){
        //Primary_TRAN
        for(int i=0;i<m.Primary_TRAN.Left.size();i++){
            m.Primary_TRAN.Left[i].erase(m.Primary_TRAN.Left[i].begin()+index);
        }
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==index){
                m.Primary_TRAN.Left.erase(m.Primary_TRAN.Left.begin()+i);
                m.Primary_TRAN.Right.erase(m.Primary_TRAN.Right.begin()+i);
                break;
            }
        }
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index>index){
                m.Primary_TRAN.Right[i].index--;
            }
        }
        //TRAN
        for(int i=0;i<m.TRAN.Left.size();i++){
            m.TRAN.Left[i].erase(m.TRAN.Left[i].begin()+index);
        }
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==index){
                m.TRAN.Left.erase(m.TRAN.Left.begin()+i);
                m.TRAN.Right.erase(m.TRAN.Right.begin()+i);
                break;
            }
        }
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index>index){
                m.TRAN.Right[i].index--;
            }
        }
        //DC
        for(int i=0;i<m.DC.Left.size();i++){
            m.DC.Left[i].erase(m.DC.Left[i].begin()+index);
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==index){
                m.DC.Left.erase(m.DC.Left.begin()+i);
                m.DC.Right.erase(m.DC.Right.begin()+i);
                break;
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index>index){
                m.DC.Right[i].index--;
            }
        }

    }
    void Delete_GND(Matrix_solve &m){
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==index){
                m.Primary_TRAN.Left[i][index]=0.0;
            }
        }
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==index){
                m.TRAN.Left[i][index]=0.0;
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==index){
                m.DC.Left[i][index]=0.0;
            }
        }
    }
};
class Element{
protected:
    string type;
    string name;
    string value;
public:
    int Right_index1=0,Right_index2=0;
    int TRAN_current_index=0;
    int Primary_current_index=0;
    int DC_current_index=0;
    Node node1,node2;
    Element(string t,string n, string v) : type(t),name(n),value(v){}
    virtual void Add_Equation(Matrix_solve &m,int node_size)=0;
    virtual void Delete_Equation(Matrix_solve &m,int node_size)=0;
    virtual void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep)=0;
    string getName(){
        return name;
    }
    void setValue(string v){
        if(stod(v)>0)value=v;
    }
    string showValue(){
        return value;
    }
    double getValue(){
        double v= stod(value);
        if(value.find("Meg")!=-1||value.find("M")!=-1)v*=1e+6;
        else if(value.find("k")!=-1)v*=1000.0;
        else if(value.find("u")!=-1)v*=1e-6;
        else if(value.find("n")!=-1)v*=1e-9;
        else if(value.find("m")!=-1)v*=0.001;
        return v;
    };
    string getType(){
        return type;
    }
};
class Resistor : public Element{
public:
    Resistor(string t,string n, string v) : Element(t,n,v){}
    void Add_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Left[i][node1.index]-=1.0/getValue();
                m.Primary_TRAN.Left[i][node2.index]+=1.0/getValue();
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Left[i][node2.index]-=1.0/getValue();
                m.Primary_TRAN.Left[i][node1.index]+=1.0/getValue();
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Left[i][node1.index]-=1.0/getValue();
                m.DC.Left[i][node2.index]+=1.0/getValue();
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node1.is_ground){
                m.DC.Left[i][node2.index]-=1.0/getValue();
                m.DC.Left[i][node1.index]+=1.0/getValue();
            }
        }
    }
    void Delete_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Left[i][node1.index]+=1.0/getValue();
                m.Primary_TRAN.Left[i][node2.index]-=1.0/getValue();
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Left[i][node2.index]+=1.0/getValue();
                m.Primary_TRAN.Left[i][node1.index]-=1.0/getValue();
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Left[i][node1.index]+=1.0/getValue();
                m.DC.Left[i][node2.index]-=1.0/getValue();
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node1.is_ground){
                m.DC.Left[i][node2.index]+=1.0/getValue();
                m.DC.Left[i][node1.index]-=1.0/getValue();
            }
        }
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Left[i][node1.index]-=1.0/getValue();
                m.TRAN.Left[i][node2.index]+=1.0/getValue();
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Left[i][node2.index]-=1.0/getValue();
                m.TRAN.Left[i][node1.index]+=1.0/getValue();
            }
        }
    }
};
class Capacitor : public Element{
public:
    double previous_current,previous_voltage;
    Capacitor(string t,string n, string v) : Element(t,n,v){}
    void Add_Equation(Matrix_solve &m,int node_size) override{
        Primary_current_index=m.Primary_TRAN.Left.size()-node_size;
        for(int i=0;i<m.Primary_TRAN.Left.size();i++){
            m.Primary_TRAN.Left[i].push_back(0.0);
        }
        vector<double> v(m.Primary_TRAN.Left.size()+1,0.0);
        m.Primary_TRAN.Left.push_back(v);
        m.Primary_TRAN.Right.push_back({0.0,0,Primary_current_index});
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][node1.index]+=1.0;
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][node2.index]-=1.0;
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Left[i][m.Primary_TRAN.Left[i].size()-1]-=1.0;
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Left[i][m.Primary_TRAN.Left[i].size()-1]+=1.0;
            }
        }
    }
    void Delete_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.Primary_TRAN.Left.size();i++){
            m.Primary_TRAN.Left[i].erase(m.Primary_TRAN.Left[i].begin()+Primary_current_index+node_size);
        }
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(!m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==Primary_current_index){
                m.Primary_TRAN.Right.erase(m.Primary_TRAN.Right.begin()+i);
                m.Primary_TRAN.Left.erase(m.Primary_TRAN.Left.begin()+i);
                break;
            }
        }
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(!m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index>Primary_current_index)m.Primary_TRAN.Right[i].index--;
        }
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Left[i][node1.index]-=2*getValue()/tstep;
                m.TRAN.Left[i][node2.index]+=2*getValue()/tstep;
                m.TRAN.Right[i].value-=(2*getValue()*previous_voltage/tstep+previous_current);
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Left[i][node1.index]+=2*getValue()/tstep;
                m.TRAN.Left[i][node2.index]-=2*getValue()/tstep;
                m.TRAN.Right[i].value+=(2*getValue()*previous_voltage/tstep+previous_current);
            }
        }
    }
};
class Inductor : public Element{
public:
    double previous_current,previous_voltage;
    Inductor(string t,string n, string v) : Element(t,n,v){}
    void Add_Equation(Matrix_solve &m,int node_size) override{
        DC_current_index=m.DC.Left.size()-node_size;
        for(int i=0;i<m.DC.Left.size();i++){
            m.DC.Left[i].push_back(0.0);
        }
        vector<double> v(m.DC.Left.size()+1,0.0);
        m.DC.Left.push_back(v);
        m.DC.Right.push_back({0.0,0,DC_current_index});
        m.DC.Left[m.DC.Left.size()-1][node1.index]+=1.0;
        m.DC.Left[m.DC.Left.size()-1][node2.index]-=1.0;
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Left[i][m.DC.Left[i].size()-1]-=1.0;
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node2.is_ground){
                m.DC.Left[i][m.DC.Left[i].size()-1]+=1.0;
            }
        }
    }
    void Delete_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.DC.Left.size();i++){
            m.DC.Left[i].push_back(0.0);
            m.DC.Left[i].erase(m.DC.Left[i].begin()+DC_current_index+node_size);
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(!m.DC.Right[i].is_node&&m.DC.Right[i].index==DC_current_index){
                m.DC.Right.erase(m.DC.Right.begin()+i);
                m.DC.Left.erase(m.DC.Left.begin()+i);
                break;
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(!m.DC.Right[i].is_node&&m.DC.Right[i].index>DC_current_index)m.DC.Right[i].index--;
        }
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Left[i][node1.index]-=tstep/(2*getValue());
                m.TRAN.Left[i][node2.index]+=tstep/(2*getValue());
                m.TRAN.Right[i].value-=(tstep*previous_voltage/(2*getValue())+previous_current);
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Left[i][node1.index]+=tstep/(2*getValue());
                m.TRAN.Left[i][node2.index]-=tstep/(2*getValue());
                m.TRAN.Right[i].value+=(tstep*previous_voltage/(2*getValue())+previous_current);
            }
        }
    }
};
class Diode : public Element{
public:
    int current_index;
    Diode(string t,string n, string v) : Element(t,n,v){}
    void Add_Equation(Matrix_solve &m,int node_size) override{}
    void Delete_Equation(Matrix_solve &m,int node_size) override{}
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{}
};
class VoltageSource : public Element{
public:
    VoltageSource(string t,string n, string v) : Element(t,n,v){}
    void Add_Equation(Matrix_solve &m,int node_size) override{
        Primary_current_index=m.Primary_TRAN.Left.size()-node_size;
        for(int i=0;i<m.Primary_TRAN.Left.size();i++){
            m.Primary_TRAN.Left[i].push_back(0.0);
        }
        vector<double> v(m.Primary_TRAN.Left.size()+1,0.0);
        m.Primary_TRAN.Left.push_back(v);
        m.Primary_TRAN.Right.push_back({getValue(),0,Primary_current_index});
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][node1.index]+=1.0;
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][node2.index]-=1.0;
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Left[i][m.Primary_TRAN.Left[i].size()-1]-=1.0;
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Left[i][m.Primary_TRAN.Left[i].size()-1]+=1.0;
            }
        }
        DC_current_index=m.DC.Left.size()-node_size;
        for(int i=0;i<m.DC.Left.size();i++){
            m.DC.Left[i].push_back(0.0);
        }
        v.resize(m.DC.Left.size()+1,0.0);
        m.DC.Left.push_back(v);
        m.DC.Right.push_back({getValue(),0,DC_current_index});
        m.DC.Left[m.DC.Left.size()-1][node1.index]+=1.0;
        m.DC.Left[m.DC.Left.size()-1][node2.index]-=1.0;
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Left[i][m.DC.Left[i].size()-1]-=1.0;
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node2.is_ground){
                m.DC.Left[i][m.DC.Left[i].size()-1]+=1.0;
            }
        }
    }
    void Delete_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.Primary_TRAN.Left.size();i++){
            m.Primary_TRAN.Left[i].erase(m.Primary_TRAN.Left[i].begin()+Primary_current_index+node_size);
        }
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(!m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==Primary_current_index){
                m.Primary_TRAN.Right.erase(m.Primary_TRAN.Right.begin()+i);
                m.Primary_TRAN.Left.erase(m.Primary_TRAN.Left.begin()+i);
                break;
            }
        }
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(!m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index>Primary_current_index)m.Primary_TRAN.Right[i].index--;
        }

        for(int i=0;i<m.DC.Left.size();i++){
            m.DC.Left[i].push_back(0.0);
            m.DC.Left[i].erase(m.DC.Left[i].begin()+DC_current_index+node_size);
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(!m.DC.Right[i].is_node&&m.DC.Right[i].index==DC_current_index){
                m.DC.Right.erase(m.DC.Right.begin()+i);
                m.DC.Left.erase(m.DC.Left.begin()+i);
                break;
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(!m.DC.Right[i].is_node&&m.DC.Right[i].index>DC_current_index)m.DC.Right[i].index--;
        }
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        TRAN_current_index=m.TRAN.Left.size()-node_size;
        for(int i=0;i<m.TRAN.Left.size();i++){
            m.TRAN.Left[i].push_back(0.0);
        }
        vector<double> v(m.TRAN.Left.size()+1,0.0);
        m.TRAN.Left.push_back(v);
        m.TRAN.Right.push_back({getValue(),0,TRAN_current_index});
        m.TRAN.Left[m.TRAN.Left.size()-1][node1.index]+=1.0;
        m.TRAN.Left[m.TRAN.Left.size()-1][node2.index]-=1.0;
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Left[i][m.TRAN.Left[i].size()-1]-=1.0;
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Left[i][m.TRAN.Left[i].size()-1]+=1.0;
            }
        }
    }
};
class CurrentSource : public Element{
public:
    CurrentSource(string t,string n, string v) : Element(t,n,v){}
    void Add_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Right[i].value+=getValue();
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Right[i].value-=getValue();
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Right[i].value+=getValue();
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node2.is_ground){
                m.DC.Right[i].value-=getValue();
            }
        }
    }
    void Delete_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Right[i].value-=getValue();
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Right[i].value+=getValue();
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Right[i].value-=getValue();
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node2.is_ground){
                m.DC.Right[i].value+=getValue();
            }
        }
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Right[i].value+=getValue();
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Right[i].value-=getValue();
            }
        }
    }
};
class Vsin : public VoltageSource{
private:
    double Vamplitude,Frequency;
public:
    Vsin(string t,string n,string v , double vamp, double freq): VoltageSource(t,n,v),Vamplitude(vamp),Frequency(freq){}
    double get_voltage(double time){
        return getValue()+Vamplitude*sin(2.0*3.1415926536*Frequency*time);
    }
    void View_profile(){
        cout << "Element name : " << getName();
        cout << " _ Nodes : " << node1.name << ", " << node2.name;
        cout << " Type : " << "SIN(" << value << ", " << Vamplitude << ", " << Frequency << ")" << endl;
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        TRAN_current_index=m.TRAN.Left.size()-node_size;
        for(int i=0;i<m.TRAN.Left.size();i++){
            m.TRAN.Left[i].push_back(0.0);
        }
        vector<double> v(m.TRAN.Left.size()+1,0.0);
        m.TRAN.Left.push_back(v);
        m.TRAN.Right.push_back({get_voltage(tstep),0,TRAN_current_index});
        m.TRAN.Left[m.TRAN.Left.size()-1][node1.index]+=1.0;
        m.TRAN.Left[m.TRAN.Left.size()-1][node2.index]-=1.0;
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Left[i][m.TRAN.Left[i].size()-1]-=1.0;
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Left[i][m.TRAN.Left[i].size()-1]+=1.0;
            }
        }
    }
};
class Isin : public CurrentSource{
private:
    double Iamplitude,Frequency;
public:
    Isin(string t,string n,string v , double iamp, double freq): CurrentSource(t,n, v),Iamplitude(iamp),Frequency(freq){}
    double get_current(double time){
        return getValue()+Iamplitude*sin(2.0*3.1415926536*Frequency*time);
    }
    void View_profile(){
        cout << "Element name : " << getName();
        cout << " _ Nodes : " << node1.name << ", " << node2.name;
        cout << " Type : " << "SIN(" << value << ", " << Iamplitude << ", " << Frequency << ")" << endl;
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Right[i].value+= get_current(tstep);
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Right[i].value-=get_current(tstep);
            }
        }
    }
};
class V_v : public VoltageSource{
public:
    Node  cntr_node1,cntr_node2;
    V_v(string t,string n, string v,Node cn1,Node cn2): VoltageSource(t,n,v),cntr_node1(cn1),cntr_node2(cn2){}
    void Add_Equation(Matrix_solve &m,int node_size) override{
        Primary_current_index=m.Primary_TRAN.Left.size()-node_size;
        for(int i=0;i<m.Primary_TRAN.Left.size();i++){
            m.Primary_TRAN.Left[i].push_back(0.0);
        }
        vector<double> v(m.Primary_TRAN.Left.size()+1,0.0);
        m.Primary_TRAN.Left.push_back(v);
        m.Primary_TRAN.Right.push_back({0.0,0,Primary_current_index});
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][node1.index]+=1.0;
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][node2.index]-=1.0;
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][cntr_node1.index]-=getValue();
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][cntr_node2.index]+=getValue();
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Left[i][m.Primary_TRAN.Left[i].size()-1]-=1.0;
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Left[i][m.Primary_TRAN.Left[i].size()-1]+=1.0;
            }
        }
        DC_current_index=m.DC.Left.size()-node_size;
        for(int i=0;i<m.DC.Left.size();i++){
            m.DC.Left[i].push_back(0.0);
        }
        v.resize(m.DC.Left.size()+1,0.0);
        m.DC.Left.push_back(v);
        m.DC.Right.push_back({0.0,0,DC_current_index});
        m.DC.Left[m.DC.Left.size()-1][node1.index]+=1.0;
        m.DC.Left[m.DC.Left.size()-1][node2.index]-=1.0;
        m.DC.Left[m.DC.Left.size()-1][cntr_node1.index]-=getValue();
        m.DC.Left[m.DC.Left.size()-1][cntr_node2.index]+=getValue();
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Left[i][m.DC.Left[i].size()-1]-=1.0;
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node2.is_ground){
                m.DC.Left[i][m.DC.Left[i].size()-1]+=1.0;
            }
        }
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        TRAN_current_index=m.TRAN.Left.size()-node_size;
        for(int i=0;i<m.TRAN.Left.size();i++){
            m.TRAN.Left[i].push_back(0.0);
        }
        vector<double> v(m.TRAN.Left.size()+1,0.0);
        m.TRAN.Left.push_back(v);
        m.TRAN.Right.push_back({getValue(),0,TRAN_current_index});
        m.TRAN.Left[m.TRAN.Left.size()-1][node1.index]+=1.0;
        m.TRAN.Left[m.TRAN.Left.size()-1][node2.index]-=1.0;
        m.TRAN.Left[m.TRAN.Left.size()-1][cntr_node1.index]-=getValue();
        m.TRAN.Left[m.TRAN.Left.size()-1][cntr_node2.index]+=getValue();
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Left[i][m.TRAN.Left[i].size()-1]-=1.0;
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Left[i][m.TRAN.Left[i].size()-1]+=1.0;
            }
        }
    }
};
class I_v : public CurrentSource{
public:
    Node cntr_node1,cntr_node2;
    I_v(string t,string n, string v,Node cn1,Node cn2): CurrentSource(t,n,v),cntr_node1(cn1),cntr_node2(cn2){}
    void Add_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Left[i][cntr_node1.index]-=getValue();
                m.Primary_TRAN.Left[i][cntr_node2.index]+=getValue();
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Left[i][cntr_node1.index]+=getValue();
                m.Primary_TRAN.Left[i][cntr_node2.index]-=getValue();
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Left[i][cntr_node1.index]-=getValue();
                m.DC.Left[i][cntr_node2.index]+=getValue();
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node2.is_ground){
                m.DC.Left[i][cntr_node1.index]+=getValue();
                m.DC.Left[i][cntr_node2.index]-=getValue();
            }
        }
    }
    void Delete_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Left[i][cntr_node1.index]+=getValue();
                m.Primary_TRAN.Left[i][cntr_node2.index]-=getValue();
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Left[i][cntr_node1.index]-=getValue();
                m.Primary_TRAN.Left[i][cntr_node2.index]+=getValue();
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Left[i][cntr_node1.index]+=getValue();
                m.DC.Left[i][cntr_node2.index]-=getValue();
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node2.is_ground){
                m.DC.Left[i][cntr_node1.index]-=getValue();
                m.DC.Left[i][cntr_node2.index]+=getValue();
            }
        }
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Left[i][cntr_node1.index]-=getValue();
                m.TRAN.Left[i][cntr_node2.index]+=getValue();
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Left[i][cntr_node1.index]+=getValue();
                m.TRAN.Left[i][cntr_node2.index]-=getValue();
            }
        }
    }
};
class V_i : public VoltageSource{
    int Cntr_TRAN_Current_index;
    int Cntr_Primary_Current_index;
    int Cntr_DC_Current_index;
public:
    string cntr_element;
    void get_cntr_current_indexes(int t,int p,int d){
        Cntr_TRAN_Current_index=t;
        Cntr_Primary_Current_index=p;
        Cntr_DC_Current_index=d;
    }
    V_i(string t,string n, string v,string ce): VoltageSource(t,n,v),cntr_element(ce){}
    void Add_Equation(Matrix_solve &m,int node_size) override{
        Primary_current_index=m.Primary_TRAN.Left.size()-node_size;
        for(int i=0;i<m.Primary_TRAN.Left.size();i++){
            m.Primary_TRAN.Left[i].push_back(0.0);
        }
        vector<double> v(m.Primary_TRAN.Left.size()+1,0.0);
        m.Primary_TRAN.Left.push_back(v);
        m.Primary_TRAN.Right.push_back({0.0,0,Primary_current_index});
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][node1.index]+=1.0;
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][node2.index]-=1.0;
        m.Primary_TRAN.Left[m.Primary_TRAN.Left.size()-1][Cntr_Primary_Current_index+node_size]-=getValue();
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Left[i][m.Primary_TRAN.Left[i].size()-1]-=1.0;
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Left[i][m.Primary_TRAN.Left[i].size()-1]+=1.0;
            }
        }
        DC_current_index=m.DC.Left.size()-node_size;
        for(int i=0;i<m.DC.Left.size();i++){
            m.DC.Left[i].push_back(0.0);
        }
        v.resize(m.DC.Left.size()+1,0.0);
        m.DC.Left.push_back(v);
        m.DC.Right.push_back({0.0,0,DC_current_index});
        m.DC.Left[m.DC.Left.size()-1][node1.index]+=1.0;
        m.DC.Left[m.DC.Left.size()-1][node2.index]-=1.0;
        m.DC.Left[m.DC.Left.size()-1][Cntr_DC_Current_index+node_size]-=getValue();
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Left[i][m.DC.Left[i].size()-1]-=1.0;
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node2.is_ground){
                m.DC.Left[i][m.DC.Left[i].size()-1]+=1.0;
            }
        }
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        TRAN_current_index=m.TRAN.Left.size()-node_size;
        for(int i=0;i<m.TRAN.Left.size();i++){
            m.TRAN.Left[i].push_back(0.0);
        }
        vector<double> v(m.TRAN.Left.size()+1,0.0);
        m.TRAN.Left.push_back(v);
        m.TRAN.Right.push_back({0.0,0,TRAN_current_index});
        m.TRAN.Left[m.TRAN.Left.size()-1][node1.index]+=1.0;
        m.TRAN.Left[m.TRAN.Left.size()-1][node2.index]-=1.0;
        m.TRAN.Left[m.TRAN.Left.size()-1][Cntr_TRAN_Current_index+node_size]-=getValue();
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Left[i][m.TRAN.Left[i].size()-1]-=1.0;
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Left[i][m.TRAN.Left[i].size()-1]+=1.0;
            }
        }
    }
};
class I_i : public CurrentSource{
    int Cntr_TRAN_Current_index;
    int Cntr_Primary_Current_index;
    int Cntr_DC_Current_index;
public:
    string cntr_element;
    void get_cntr_current_indexes(int t,int p,int d){
        Cntr_TRAN_Current_index=t;
        Cntr_Primary_Current_index=p;
        Cntr_DC_Current_index=d;
    }
    I_i(string t,string n, string v,string ce): CurrentSource(t,n,v),cntr_element(ce){}
    void Add_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Left[i][Cntr_Primary_Current_index+node_size]-=getValue();
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Left[i][Cntr_Primary_Current_index+node_size]+=getValue();
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Left[i][Cntr_DC_Current_index+node_size]-=getValue();
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node2.is_ground){
                m.DC.Left[i][Cntr_DC_Current_index+node_size]+=getValue();
            }
        }
    }
    void Delete_Equation(Matrix_solve &m,int node_size) override{
        for(int i=0;i<m.Primary_TRAN.Right.size();i++){
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.Primary_TRAN.Left[i][Cntr_Primary_Current_index+node_size]+=getValue();
            }
            if(m.Primary_TRAN.Right[i].is_node&&m.Primary_TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.Primary_TRAN.Left[i][Cntr_Primary_Current_index+node_size]-=getValue();
            }
        }
        for(int i=0;i<m.DC.Right.size();i++){
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node1.index&&!node1.is_ground){
                m.DC.Left[i][Cntr_DC_Current_index+node_size]+=getValue();
            }
            if(m.DC.Right[i].is_node&&m.DC.Right[i].index==node2.index&&!node2.is_ground){
                m.DC.Left[i][Cntr_DC_Current_index+node_size]-=getValue();
            }
        }
    }
    void Add_TRAN_Equation(Matrix_solve &m,int node_size,double tstep) override{
        for(int i=0;i<m.TRAN.Right.size();i++){
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node1.index&&!node1.is_ground){
                m.TRAN.Left[i][Cntr_TRAN_Current_index+node_size]-=getValue();
            }
            if(m.TRAN.Right[i].is_node&&m.TRAN.Right[i].index==node2.index&&!node2.is_ground){
                m.TRAN.Left[i][Cntr_TRAN_Current_index+node_size]+=getValue();
            }
        }
    }
};
class Circuit{
    bool has_changed_DC_matrix=1;
public:
    Matrix_solve matrixSolve;
    vector<unique_ptr<Element>> element={};
    vector<Node> node={};
    void containsElementWithName(bool &b, string s){
        for(int i=0;i<element.size();i++){
            if(element[i]->getName()==s){
                b=!b;
                break;
            }
        }
    }
    bool isCircuitComplete(int index){
        static vector<bool> v;
        static vector<vector<int>> NeighbourNodes={};
        if(index==0){
            NeighbourNodes.clear();
            v.clear();
            for(int i=0;i<node.size();i++) {
                NeighbourNodes.push_back({});
                v.push_back(false);
            }
            for(int i=0;i<element.size();i++){
                NeighbourNodes[element[i]->node1.index].push_back(element[i]->node2.index);
                NeighbourNodes[element[i]->node2.index].push_back(element[i]->node1.index);
            }
            for(int i=0;i<NeighbourNodes.size();i++){
                if(NeighbourNodes[i].size()==0)return false;
            }
        }
        v[index]=true;
        for(int i=0;i<NeighbourNodes[index].size();i++){
            int n=NeighbourNodes[index][i];
            if (!v[n]) isCircuitComplete(n);
        }
        if(index==0){
            for(int i=0;i<v.size();i++){
                if(!v[i]) return false;
            }
            return true;
        }
        return false;
    }
    double element_current_shower(int index,vector<double> answer,double tspent,double tstep,string analysis_type){
        string type=element[index]->getType();
        if(type=="R"){
            double v1=answer[element[index]->node1.index];
            double v2=answer[element[index]->node2.index];
            return (v1-v2)/element[index]->getValue();
        }
        else if(type=="C"){
            if(element[index]->node1.index==element[index]->node2.index)return 0.0;
            if(analysis_type=="DC")return 0.0;
            else {
                Capacitor* p=dynamic_cast<Capacitor*>(element[index].get());
                if(tspent!=0) return 2*p->getValue()*(answer[p->node1.index]-answer[p->node2.index]-p->previous_voltage)/tstep-p->previous_current;
                else return answer[p->Primary_current_index+node.size()];
            }
        }
        else if(type=="L"){
            if(element[index]->node1.index==element[index]->node2.index)return 0.0;
            Inductor* p=dynamic_cast<Inductor*>(element[index].get());
            if(analysis_type=="DC") return answer[p->DC_current_index+node.size()];
            else {
                if(tspent!=0)return tstep*(answer[p->node1.index]-answer[p->node2.index]+p->previous_voltage)/(2*p->getValue())+p->previous_current;
                else return 0.0;
            }
        }
        else if(type=="D"){
            Diode* p=dynamic_cast<Diode*>(element[index].get());
            return answer[p->current_index+node.size()];
        }
        else if(type=="V"){
            VoltageSource* p=dynamic_cast<VoltageSource*>(element[index].get());
            if(analysis_type=="DC") return answer[p->DC_current_index+node.size()];
            else {
                if(tspent==0)return answer[p->Primary_current_index+node.size()];
                return answer[p->TRAN_current_index+node.size()];
            }
        }
        else if(type=="I"){
            return element[index]->getValue();
        }
        else if(type=="Vsin"){
            Vsin* p=dynamic_cast<Vsin*>(element[index].get());
            if(analysis_type=="DC") return answer[p->DC_current_index+node.size()];
            else {
                if(tspent==0)return answer[p->Primary_current_index+node.size()];
                return answer[p->TRAN_current_index+node.size()];
            }
        }
        else if(type=="Isin"){
            Isin* p = dynamic_cast<Isin*>(element[index].get());
            return p->get_current(tspent);
        }
        else if(type=="E"){
            V_v* p = dynamic_cast<V_v*>(element[index].get());
            if(analysis_type=="DC") return answer[p->DC_current_index+node.size()];
            else {
                if(tspent==0)return answer[p->Primary_current_index+node.size()];
                return answer[p->TRAN_current_index+node.size()];
            }
        }
        else if(type=="G"){
            I_v* p = dynamic_cast<I_v*>(element[index].get());
            return p->getValue()*(answer[p->cntr_node1.index]-answer[p->cntr_node2.index]);
        }
        else if(type=="H"){
            V_i* p = dynamic_cast<V_i*>(element[index].get());
            if(analysis_type=="DC") return answer[p->DC_current_index+node.size()];
            else {
                if(tspent==0)return answer[p->Primary_current_index+node.size()];
                return answer[p->TRAN_current_index+node.size()];
            }
        }
        else if(type=="F"){
            I_i* p = dynamic_cast<I_i*>(element[index].get());
            int cntr_index=0;
            for(int i=0;i<element.size();i++){
                if(p->cntr_element==element[i]->getName()){
                    cntr_index=i;
                    break;
                }
            }
            return p->getValue()* element_current_shower(cntr_index,answer,tspent,tstep,analysis_type);
        }
    }
    void Add(smatch match,string input_type){
        has_changed_DC_matrix=1;
        if(input_type=="add_GND"){
            int node_index=-1;
            for(int i=0;i<node.size();i++){
                if(node[i].name==match[2])node_index=node[i].index;
            }
            if(node_index==-1){
                node_index=node.size();
                node.push_back({match[2],node_index,1});
                node[node_index].Add_Equation(matrixSolve);
            }
            node[node_index].is_ground=1;
            node[node_index].Add_GND(matrixSolve);
            for(int i=0;i<element.size();i++){
                if(element[i]->node1.name==node[node_index].name)element[i]->node1=node[node_index];
                if(element[i]->node2.name==node[node_index].name)element[i]->node2=node[node_index];
            }
        }
        else {
            int node_index1=-1,node_index2=-1;
            double value;
            if(input_type!="add_Diode"){
                for(int i=0;i<node.size();i++){
                    if(node[i].name==match[2])node_index1=node[i].index;
                    if(node[i].name==match[3])node_index2=node[i].index;
                }
                if(node_index1==-1){
                    node_index1=node.size();
                    node.push_back({match[2],node_index1,0});
                    node[node_index1].Add_Equation(matrixSolve);
                }
                if(node_index2==-1){
                    node_index2=node.size();
                    node.push_back({match[3],node_index2,0});
                    node[node_index1].Add_Equation(matrixSolve);
                }
            }
            if(input_type=="add_Resistor"){
                element.push_back(make_unique<Resistor>("R",match[1],match[4]));
            }
            else if(input_type=="add_Capacitor"){
                element.push_back(make_unique<Capacitor>("C",match[1],match[4]));
            }
            else if(input_type=="add_Inductor"){
                element.push_back(make_unique<Inductor>("L",match[1],match[4]));
            }
            else if(input_type=="add_Diode"){
                if(match[4]=="D") element.push_back(make_unique<Diode>("D",match[1],"0"));
                else element.push_back(make_unique<Diode>("D",match[1],"1"));
            }
            else if(input_type=="add_VoltageSource"){
                element.push_back(make_unique<VoltageSource>("V",match[1],match[4]));
            }
            else if(input_type=="add_CurrentSource"){
                element.push_back(make_unique<CurrentSource>("I",match[1],match[4]));
            }
            else if(input_type=="add_VoltageSource_sin"){
                element.push_back(make_unique<Vsin>("Vsin",match[1],match[4], stod(match[5]), stod(match[6])));
            }
            else if(input_type=="add_CurrentSource_sin"){
                element.push_back(make_unique<Isin>("Isin",match[1],match[4], stod(match[5]), stod(match[6])));
            }
            else if(input_type=="add_VoltageSource->voltage"){
                int cntr_node1=0,cntr_node2=0;
                for(int i=0;i<node.size();i++){
                    if(node[i].name==match[4])cntr_node1=node[i].index;
                    if(node[i].name==match[5])cntr_node2=node[i].index;
                }
                element.push_back(make_unique<V_v>("E",match[1],match[6],node[cntr_node1],node[cntr_node2]));
            }
            else if(input_type=="add_CurrentSource->voltage"){
                int cntr_node1=0,cntr_node2=0;
                for(int i=0;i<node.size();i++){
                    if(node[i].name==match[4])cntr_node1=node[i].index;
                    if(node[i].name==match[5])cntr_node2=node[i].index;
                }
                element.push_back(make_unique<I_v>("G",match[1],match[6],node[cntr_node1],node[cntr_node2]));
            }
            else if(input_type=="add_VoltageSource->current"){
                element.push_back(make_unique<V_i>("H",match[1],match[5],match[4]));
                V_i* p = dynamic_cast<V_i*>(element[element.size()-1].get());
                for(int i=0;i<element.size();i++){
                    if(element[i]->getName()==match[4]){
                        if(VoltageSource* p2=dynamic_cast<VoltageSource*>(element[i].get())){
                            p->get_cntr_current_indexes(p2->TRAN_current_index,p2->Primary_current_index,p2->DC_current_index);
                            break;
                        }
                    }
                }
            }
            else if(input_type=="add_CurrentSource->current"){
                element.push_back(make_unique<I_i>("F",match[1],match[5],match[4]));
                I_i* p = dynamic_cast<I_i*>(element[element.size()-1].get());
                for(int i=0;i<element.size();i++){
                    if(element[i]->getName()==match[4]){
                        if(VoltageSource* p2=dynamic_cast<VoltageSource*>(element[i].get())){
                            p->get_cntr_current_indexes(p2->TRAN_current_index,p2->Primary_current_index,p2->DC_current_index);
                            break;
                        }
                    }
                }
            }
            if(input_type!="add_Diode") {
                element[element.size() - 1]->Add_Equation(matrixSolve, node.size());
                element[element.size() - 1]->node1 = node[node_index1];
                element[element.size() - 1]->node2 = node[node_index2];
            }
        }
    }
    void Delete(smatch match,string input_type){
        has_changed_DC_matrix=1;
        if(input_type=="delete_Element"){
            int node_index1=0,node_index2=0;
            for(int i=0;i<element.size();i++){
                if(element[i]->getName()==match[1]){
                    node_index1=element[i]->node1.index;
                    node_index2=element[i]->node2.index;
                    element[i]->Delete_Equation(matrixSolve,node.size());
                    for(int j=0;j<element.size();j++){
                        if(j!=i){
                            if(element[i]->getType()=="C"){
                                if(element[j]->Primary_current_index>element[i]->Primary_current_index)element[j]->Primary_current_index--;
                            }
                            else if(element[i]->getType()=="L"){
                                if(element[j]->DC_current_index>element[i]->DC_current_index)element[j]->DC_current_index--;
                            }
                            else if(element[i]->getType()=="V"||element[i]->getType()=="Vsin"||element[i]->getType()=="E"||element[i]->getType()=="H"){
                                if(element[j]->Primary_current_index>element[i]->Primary_current_index)element[j]->Primary_current_index--;
                                if(element[j]->DC_current_index>element[i]->DC_current_index)element[j]->DC_current_index--;
                                if(element[j]->TRAN_current_index>element[i]->TRAN_current_index)element[j]->TRAN_current_index--;
                            }
                        }
                    }
                    element.erase(element.begin()+i);
                    for(int j=0;j<element.size();j++){
                        if(V_i* p = dynamic_cast<V_i*>(element[j].get())){
                            for(int k=0;k<element.size();k++){
                                if(p->cntr_element==element[k]->getName()){
                                    p->get_cntr_current_indexes(element[k]->TRAN_current_index,element[k]->Primary_current_index,element[k]->DC_current_index);
                                    break;
                                }
                            }
                        }
                        if(I_i* p = dynamic_cast<I_i*>(element[j].get())){
                            for(int k=0;k<element.size();k++){
                                if(p->cntr_element==element[k]->getName()){
                                    p->get_cntr_current_indexes(element[k]->TRAN_current_index,element[k]->Primary_current_index,element[k]->DC_current_index);
                                    break;
                                }
                            }
                        }
                    }
                    bool is_excess_node=1,is_excess_node2=1;
                    for(int j=0;j<element.size();j++){
                        if(element[j]->node1.index==node_index1||element[j]->node2.index==node_index1) is_excess_node=0;
                        if(element[j]->node1.index==node_index2||element[j]->node2.index==node_index2) is_excess_node2=0;
                    }
                    if(is_excess_node){
                        node[node_index1].Delete_Equation(matrixSolve);
                        node.erase(node.begin()+node_index1);
                        for(int j=0;j<node.size();j++){
                            node[j].index=j;
                            for(int k=0;k<element.size();k++){
                                if(element[k]->node1.name==node[j].name)element[k]->node1=node[j];
                                if(element[k]->node2.name==node[j].name)element[k]->node2=node[j];
                            }
                        }
                        if(node_index2>node_index1)node_index2--;
                    }
                    if(is_excess_node2){
                        node[node_index2].Delete_Equation(matrixSolve);
                        node.erase(node.begin()+node_index2);
                        for(int j=0;j<node.size();j++){
                            node[j].index=j;
                            for(int k=0;k<element.size();k++){
                                if(element[k]->node1.name==node[j].name)element[k]->node1=node[j];
                                if(element[k]->node2.name==node[j].name)element[k]->node2=node[j];
                            }
                        }
                    }
                    break;
                }
            }
        }
        else if(input_type=="delete_GND"){
            for(int i=0;i<node.size();i++){
                if(node[i].name==match[1]){
                    node[i].is_ground=0;
                    node[i].Delete_GND(matrixSolve);
                    bool is_excess_node=1;
                    for(int j=0;j<element.size();j++){
                        if(element[j]->node1.name==match[1]||element[j]->node2.name==match[1]){
                            is_excess_node=0;
                            element[j]->Add_Equation(matrixSolve,node.size());
                        }
                    }
                    if(is_excess_node){
                        node[i].Delete_Equation(matrixSolve);
                        node.erase(node.begin()+i);
                        for(int j=0;j<node.size();j++){
                            node[j].index=j;
                            for(int k=0;k<element.size();k++){
                                if(element[k]->node1.name==node[j].name)element[k]->node1=node[j];
                                if(element[k]->node2.name==node[j].name)element[k]->node2=node[j];
                            }
                        }
                    }
                    break;
                }
            }
        }
    }
    void Rename(smatch match){
        for(int i=0;i<node.size();i++){
            if(node[i].name==match[1]){
                for(int j=0;j<element.size();j++){
                    if(element[j]->node1.name==match[1])element[j]->node1.name=match[2];
                    if(element[j]->node2.name==match[1])element[j]->node2.name=match[2];
                    if(V_v* p=dynamic_cast<V_v*>(element[j].get())){
                        if(p->cntr_node1.name==match[1])p->cntr_node1.name=match[2];
                        if(p->cntr_node2.name==match[1])p->cntr_node2.name=match[2];
                    }
                    if(I_v* p=dynamic_cast<I_v*>(element[j].get())){
                        if(p->cntr_node1.name==match[1])p->cntr_node1.name=match[2];
                        if(p->cntr_node2.name==match[1])p->cntr_node2.name=match[2];
                    }
                }
                node[i].name=match[2];
                break;
            }
        }
    }
    bool Print_TRAN(smatch match,double tspent){
        double tstep,tstart,tstop;
        tstep= stod(match[3]);
        tstop= stod(match[2]);
        string s_tstart=match[1];
        if(s_tstart.length()==0)s_tstart.push_back('0');
        tstart= stod(s_tstart);
        if(tstart!=0&&tspent==tstart){
            if(matrixSolve.LUsetter(matrixSolve.Primary_TRAN)){
                matrixSolve.Solve(matrixSolve.Primary_TRAN);
                matrixSolve.TRAN.Left.clear();
                matrixSolve.TRAN.Right.clear();
                vector<double> v(node.size(),0.0);
                matrixSolve.TRAN.Left.resize(node.size(),v);
                matrixSolve.TRAN.Right.resize(node.size());
                for(int i=0;i<node.size();i++){
                    matrixSolve.TRAN.Right[i].is_node= true;
                    matrixSolve.TRAN.Right[i].value=0.0;
                    matrixSolve.TRAN.Right[i].index=i;
                    if(node[i].is_ground)matrixSolve.TRAN.Left[i][i]+=1.0;
                }
                for(int i=0;i<element.size();i++){
                    if(Capacitor* p=dynamic_cast<Capacitor*>(element[i].get())){
                        p->previous_current= element_current_shower(i,matrixSolve.Primary_TRAN.Answer,0,tstep,"TRAN");
                        p->previous_voltage=matrixSolve.Primary_TRAN.Answer[p->node1.index]-matrixSolve.Primary_TRAN.Answer[p->node2.index];
                    }
                    else if(Inductor* p=dynamic_cast<Inductor*>(element[i].get())){
                        p->previous_current=0;
                        p->previous_voltage=matrixSolve.Primary_TRAN.Answer[p->node1.index]-matrixSolve.Primary_TRAN.Answer[p->node2.index];
                    }
                    element[i]->Add_TRAN_Equation(matrixSolve,node.size(),tstep);
                }
                if(matrixSolve.LUsetter(matrixSolve.TRAN)){
                    for(int i=0;i<element.size();i++){
                        if(Capacitor* p=dynamic_cast<Capacitor*>(element[i].get())){
                            for(int j=0;j<matrixSolve.TRAN.Right.size();j++){
                                if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node1.index)p->Right_index1=j;
                                if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node2.index)p->Right_index2=j;
                            }
                        }
                        else if(Inductor* p=dynamic_cast<Inductor*>(element[i].get())){
                            for(int j=0;j<matrixSolve.TRAN.Right.size();j++){
                                if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node1.index)p->Right_index1=j;
                                if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node2.index)p->Right_index2=j;
                            }
                        }
                        else if(Vsin* p=dynamic_cast<Vsin*>(element[i].get())){
                            for(int j=0;j<matrixSolve.TRAN.Right.size();j++){
                                if(!matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->TRAN_current_index)p->Right_index1=j;
                            }
                        }
                        else if(Isin* p=dynamic_cast<Isin*>(element[i].get())){
                            for(int j=0;j<matrixSolve.TRAN.Right.size();j++){
                                if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node1.index)p->Right_index1=j;
                                if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node2.index)p->Right_index2=j;
                            }
                        }
                    }
                }
                else return false;
            }
            else return false;
            int number=tstart/tstep;
            for(int j=1;j<=number;j++){
                double tspent2=double(j)*tstep;
                matrixSolve.Solve(matrixSolve.TRAN);
                for(int i=0;i<element.size();i++){
                    if(Capacitor* p=dynamic_cast<Capacitor*>(element[i].get())){
                        matrixSolve.TRAN.Right[p->Right_index1].value+=(2*p->getValue()*p->previous_voltage/tstep+p->previous_current);
                        matrixSolve.TRAN.Right[p->Right_index2].value-=(2*p->getValue()*p->previous_voltage/tstep+p->previous_current);
                        p->previous_voltage=matrixSolve.TRAN.Answer[p->node1.index]-matrixSolve.TRAN.Answer[p->node2.index];
                        p->previous_current= element_current_shower(i,matrixSolve.TRAN.Answer,tspent2,tstep,"TRAN");
                        matrixSolve.TRAN.Right[p->Right_index1].value-=(2*p->getValue()*p->previous_voltage/tstep+p->previous_current);
                        matrixSolve.TRAN.Right[p->Right_index2].value+=(2*p->getValue()*p->previous_voltage/tstep+p->previous_current);
                    }
                    else if(Inductor* p=dynamic_cast<Inductor*>(element[i].get())){
                        matrixSolve.TRAN.Right[p->Right_index1].value+=(tstep*p->previous_voltage/(2*p->getValue())+p->previous_current);
                        matrixSolve.TRAN.Right[p->Right_index2].value-=(tstep*p->previous_voltage/(2*p->getValue())+p->previous_current);
                        p->previous_voltage=matrixSolve.TRAN.Answer[p->node1.index]-matrixSolve.TRAN.Answer[p->node2.index];
                        p->previous_current= element_current_shower(i,matrixSolve.TRAN.Answer,tspent2,tstep,"TRAN");
                        matrixSolve.TRAN.Right[p->Right_index1].value-=(tstep*p->previous_voltage/(2*p->getValue())+p->previous_current);
                        matrixSolve.TRAN.Right[p->Right_index2].value+=(tstep*p->previous_voltage/(2*p->getValue())+p->previous_current);
                    }
                    else if(Vsin* p=dynamic_cast<Vsin*>(element[i].get())){
                        matrixSolve.TRAN.Right[p->Right_index1].value=p->get_voltage(tspent2);
                    }
                    else if(Isin* p=dynamic_cast<Isin*>(element[i].get())){
                        matrixSolve.TRAN.Right[p->Right_index1].value-=p->get_current(tspent2);
                        matrixSolve.TRAN.Right[p->Right_index1].value+=p->get_current(tspent2+tstep);
                        matrixSolve.TRAN.Right[p->Right_index2].value+=p->get_current(tspent2);
                        matrixSolve.TRAN.Right[p->Right_index2].value-=p->get_current(tspent2+tstep);
                    }
                }
            }
            return true;

        }
        else{
            if(tspent==0){
                if(matrixSolve.LUsetter(matrixSolve.Primary_TRAN)){
                    matrixSolve.Solve(matrixSolve.Primary_TRAN);
                    matrixSolve.TRAN.Left.clear();
                    matrixSolve.TRAN.Right.clear();
                    vector<double> v(node.size(),0.0);
                    matrixSolve.TRAN.Left.resize(node.size(),v);
                    matrixSolve.TRAN.Right.resize(node.size());
                    for(int i=0;i<node.size();i++){
                        matrixSolve.TRAN.Right[i].is_node= true;
                        matrixSolve.TRAN.Right[i].value=0.0;
                        matrixSolve.TRAN.Right[i].index=i;
                        if(node[i].is_ground)matrixSolve.TRAN.Left[i][i]+=1.0;
                    }
                    for(int i=0;i<element.size();i++){
                        if(Capacitor* p=dynamic_cast<Capacitor*>(element[i].get())){
                            p->previous_current= element_current_shower(i,matrixSolve.Primary_TRAN.Answer,0,tstep,"TRAN");
                            p->previous_voltage=matrixSolve.Primary_TRAN.Answer[p->node1.index]-matrixSolve.Primary_TRAN.Answer[p->node2.index];
                        }
                        else if(Inductor* p=dynamic_cast<Inductor*>(element[i].get())){
                            p->previous_current=0;
                            p->previous_voltage=matrixSolve.Primary_TRAN.Answer[p->node1.index]-matrixSolve.Primary_TRAN.Answer[p->node2.index];
                        }
                        element[i]->Add_TRAN_Equation(matrixSolve,node.size(),tstep);
                    }
                    if(matrixSolve.LUsetter(matrixSolve.TRAN)){
                        for(int i=0;i<element.size();i++){
                            if(Capacitor* p=dynamic_cast<Capacitor*>(element[i].get())){
                                for(int j=0;j<matrixSolve.TRAN.Right.size();j++){
                                    if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node1.index)p->Right_index1=j;
                                    if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node2.index)p->Right_index2=j;
                                }
                            }
                            else if(Inductor* p=dynamic_cast<Inductor*>(element[i].get())){
                                for(int j=0;j<matrixSolve.TRAN.Right.size();j++){
                                    if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node1.index)p->Right_index1=j;
                                    if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node2.index)p->Right_index2=j;
                                }
                            }
                            else if(Vsin* p=dynamic_cast<Vsin*>(element[i].get())){
                                for(int j=0;j<matrixSolve.TRAN.Right.size();j++){
                                    if(!matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->TRAN_current_index)p->Right_index1=j;
                                }
                            }
                            else if(Isin* p=dynamic_cast<Isin*>(element[i].get())){
                                for(int j=0;j<matrixSolve.TRAN.Right.size();j++){
                                    if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node1.index)p->Right_index1=j;
                                    if(matrixSolve.TRAN.Right[j].is_node&&matrixSolve.TRAN.Right[j].index==p->node2.index)p->Right_index2=j;
                                }
                            }
                        }
                        return true;
                    }
                    else return false;
                }
                else return false;
            }
            else {
                matrixSolve.Solve(matrixSolve.TRAN);
                for(int i=0;i<element.size();i++){
                    if(Capacitor* p=dynamic_cast<Capacitor*>(element[i].get())){
                        matrixSolve.TRAN.Right[p->Right_index1].value+=(2*p->getValue()*p->previous_voltage/tstep+p->previous_current);
                        matrixSolve.TRAN.Right[p->Right_index2].value-=(2*p->getValue()*p->previous_voltage/tstep+p->previous_current);
                        p->previous_voltage=matrixSolve.TRAN.Answer[p->node1.index]-matrixSolve.TRAN.Answer[p->node2.index];
                        p->previous_current= element_current_shower(i,matrixSolve.TRAN.Answer,tspent,tstep,"TRAN");
                        matrixSolve.TRAN.Right[p->Right_index1].value-=(2*p->getValue()*p->previous_voltage/tstep+p->previous_current);
                        matrixSolve.TRAN.Right[p->Right_index2].value+=(2*p->getValue()*p->previous_voltage/tstep+p->previous_current);
                    }
                    else if(Inductor* p=dynamic_cast<Inductor*>(element[i].get())){
                        matrixSolve.TRAN.Right[p->Right_index1].value+=(tstep*p->previous_voltage/(2*p->getValue())+p->previous_current);
                        matrixSolve.TRAN.Right[p->Right_index2].value-=(tstep*p->previous_voltage/(2*p->getValue())+p->previous_current);
                        p->previous_voltage=matrixSolve.TRAN.Answer[p->node1.index]-matrixSolve.TRAN.Answer[p->node2.index];
                        p->previous_current= element_current_shower(i,matrixSolve.TRAN.Answer,tspent,tstep,"TRAN");
                        matrixSolve.TRAN.Right[p->Right_index1].value-=(tstep*p->previous_voltage/(2*p->getValue())+p->previous_current);
                        matrixSolve.TRAN.Right[p->Right_index2].value+=(tstep*p->previous_voltage/(2*p->getValue())+p->previous_current);
                    }
                    else if(Vsin* p=dynamic_cast<Vsin*>(element[i].get())){
                        matrixSolve.TRAN.Right[p->Right_index1].value=p->get_voltage(tspent);
                    }
                    else if(Isin* p=dynamic_cast<Isin*>(element[i].get())){
                        matrixSolve.TRAN.Right[p->Right_index1].value-=p->get_current(tspent);
                        matrixSolve.TRAN.Right[p->Right_index1].value+=p->get_current(tspent+tstep);
                        matrixSolve.TRAN.Right[p->Right_index2].value+=p->get_current(tspent);
                        matrixSolve.TRAN.Right[p->Right_index2].value-=p->get_current(tspent+tstep);
                    }
                }
                return true;
            }
        }
    }
    bool Print_DC(double final_value,int source_index){
        if(has_changed_DC_matrix){
            if(matrixSolve.LUsetter(matrixSolve.DC)) has_changed_DC_matrix=0;
            else return false;
        }
        double first_value=element[source_index]->getValue();
        if(element[source_index]->getType()=="V"){
            for(int i=0;i<matrixSolve.DC.Right.size();i++){
                if(!matrixSolve.DC.Right[i].is_node&&matrixSolve.DC.Right[i].index==element[source_index]->DC_current_index){
                    matrixSolve.DC.Right[i].value-=first_value;
                    matrixSolve.DC.Right[i].value+=final_value;
                    matrixSolve.Solve(matrixSolve.DC);
                    matrixSolve.DC.Right[i].value+=first_value;
                    matrixSolve.DC.Right[i].value-=final_value;
                    return true;
                }
            }
        }
        else if(element[source_index]->getType()=="I"){
            for(int i=0;i<matrixSolve.DC.Right.size();i++){
                if(matrixSolve.DC.Right[i].is_node&&matrixSolve.DC.Right[i].index==element[source_index]->node1.index&&!element[source_index]->node1.is_ground){
                    matrixSolve.DC.Right[i].value-=first_value;
                    matrixSolve.DC.Right[i].value+=final_value;
                }
                if(matrixSolve.DC.Right[i].is_node&&matrixSolve.DC.Right[i].index==element[source_index]->node2.index&&!element[source_index]->node2.is_ground){
                    matrixSolve.DC.Right[i].value+=first_value;
                    matrixSolve.DC.Right[i].value-=final_value;
                }
            }
            matrixSolve.Solve(matrixSolve.DC);
            for(int i=0;i<matrixSolve.DC.Right.size();i++){
                if(matrixSolve.DC.Right[i].is_node&&matrixSolve.DC.Right[i].index==element[source_index]->node1.index&&!element[source_index]->node1.is_ground){
                    matrixSolve.DC.Right[i].value+=first_value;
                    matrixSolve.DC.Right[i].value-=final_value;
                }
                if(matrixSolve.DC.Right[i].is_node&&matrixSolve.DC.Right[i].index==element[source_index]->node2.index&&!element[source_index]->node2.is_ground){
                    matrixSolve.DC.Right[i].value-=first_value;
                    matrixSolve.DC.Right[i].value+=final_value;
                }
            }
            return true;
        }
    }
};
class View{
private:
    Circuit circuit;
    string input_type="";
    smatch match;
    string get="";
    void Error_handling(){
        if(regex_match(get,regex (R"(^\s*add\b.*)"))) {
            regex pattern[]={regex(R"(^\s*add\s+([^CLDVI]\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*$)"),
                             regex(R"(^\s*add\s+([^RLDVI]\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*$)"),
                             regex(R"(^\s*add\s+([^RCDVI]\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*$)"),
                             regex(R"(^\s*add\s+([^RCLVI]\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s*$)"),
                             regex(R"(^\s*add\s+(\w+)\s+(\w+)\s*$)"),
                             regex(R"(^\s*add\s+([^RCLDI]\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*$)"),
                             regex(R"(^\s*add\s+([^RCLDV]\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*$)"),
                             regex(R"(^\s*add\s+([^I]\w+)\s+(\w+)\s+(\w+)\s+SIN\(\s*([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s+([-+]?\d+(?:\.\d+)?(?:e[-+]?\d+(?:\.\d+)?)?)\s+([-+]?\d+(?:\.\d+)?(?:e[-+]?\d+(?:\.\d+)?)?)\s*\)\s*$)"),
                             regex(R"(^\s*add\s+([^V]\w+)\s+(\w+)\s+(\w+)\s+SIN\(\s*([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s+([-+]?\d+(?:\.\d+)?(?:e[-+]?\d+(?:\.\d+)?)?)\s+([-+]?\d+(?:\.\d+)?(?:e[-+]?\d+(?:\.\d+)?)?)\s*\)\s*$)"),
                             regex(R"(^\s*add\s+([^G]\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)"),//E
                             regex(R"(^\s*add\s+([^E]\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)"),//G
                             regex(R"(^\s*add\s+([^F]\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)"),//H
                             regex(R"(^\s*add\s+([^H]\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)")//F
            };
            if(regex_search(get,match,pattern[0])){
                if(match[1].str()[0]!='R') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else if(stod(match[4])<=0) cout << "Error: Resistance cannot be zero or negative" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: Resistor " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_Resistor";
                }
            }
            else if(regex_search(get,match,pattern[1])){
                if(match[1].str()[0]!='C') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else if(stod(match[4])<=0) cout << "Error: Capacitance cannot be zero or negative" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: Capacitor " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_Capacitor";
                }
            }
            else if(regex_search(get,match,pattern[2])){
                if(match[1].str()[0]!='L') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else if(stod(match[4])<=0) cout << "Error: Inductance cannot be zero or negative" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: inductor " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_Inductor";
                }
            }
            else if(regex_search(get,match,pattern[3])){
                if(match[1].str()[0]!='D') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else if(match[4]!="D"&&match[4]!="Z") cout << "Error: Model " << match[4] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: diode " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_Diode";
                }
            }
            else if(regex_search(get,match,pattern[4])){
                if(match[1]!="GND") cout << "Error: Element GND not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else{
                    bool first_ground=1;
                    for(int i=0;i<circuit.node.size();i++){
                        if(circuit.node[i].is_ground){
                            first_ground=0;
                            break;
                        }
                    }
                    if(!first_ground) cout << "Error: GND already added" << endl;
                    else input_type="add_GND";
                }
            }
            else if(regex_search(get,match,pattern[5])){
                if(match[1].str()[0]!='V') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: VoltageSource " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_VoltageSource";
                }
            }
            else if(regex_search(get,match,pattern[6])){
                if(match[1].str()[0]!='I') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: CurrentSource " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_CurrentSource";
                }
            }
            else if(regex_search(get,match,pattern[7])){
                if(match[1].str()[0]!='V') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: VoltageSource " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_VoltageSource_sin";
                }
            }
            else if(regex_search(get,match,pattern[8])){
                if(match[1].str()[0]!='I') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: CurrentSource " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_CurrentSource_sin";
                }
            }
            else if(regex_search(get,match,pattern[9])){
                if(match[1].str()[0]!='E') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else if(match[4].str()[0]!='N'||match[4].length()==1)cout << "Error: Node " << match[4] << " not found in library" << endl;
                else if(match[5].str()[0]!='N'||match[5].length()==1)cout << "Error: Node " << match[5] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    int exist_control=0;
                    circuit.containsElementWithName(new_element,match[1]);
                    for(int i=0;i<circuit.node.size();i++){
                        if(circuit.node[i].name==match[4]||circuit.node[i].name==match[5]){
                            if(exist_control==0)exist_control=1;
                            else if(exist_control==1){
                                exist_control=2;
                                break;
                            }
                        }
                    }
                    if(!new_element) cout << "Error: VoltageSource " << match[1] << " already exists in the circuit" << endl;
                    else if(exist_control!=2) cout  << "Error: Dependent source has an undefined control element." << endl;
                    else input_type="add_VoltageSource->voltage";
                }
            }
            else if(regex_search(get,match,pattern[10])){
                if(match[1].str()[0]!='G') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else if(match[4].str()[0]!='N'||match[4].length()==1)cout << "Error: Node " << match[4] << " not found in library" << endl;
                else if(match[5].str()[0]!='N'||match[5].length()==1)cout << "Error: Node " << match[5] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    int exist_control=0;
                    circuit.containsElementWithName(new_element,match[1]);
                    for(int i=0;i<circuit.node.size();i++){
                        if(circuit.node[i].name==match[4]||circuit.node[i].name==match[5]){
                            if(exist_control==0)exist_control=1;
                            else if(exist_control==1){
                                exist_control=2;
                                break;
                            }
                        }
                    }
                    if(!new_element) cout << "Error: CurrentSource " << match[1] << " already exists in the circuit" << endl;
                    else if(exist_control!=2) cout  << "Error: Dependent source has an undefined control element." << endl;
                    else input_type="add_CurrentSource->voltage";
                }
            }
            else if(regex_search(get,match,pattern[11])){
                if(match[1].str()[0]!='H') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else {
                    bool new_element=1,exist_control=0;
                    circuit.containsElementWithName(new_element,match[1]);
                    for(int i=0;i<circuit.element.size();i++){
                        if(circuit.element[i]->getType()=="V"&&circuit.element[i]->getName()==match[4]){
                            exist_control=1;
                            break;
                        }
                    }
                    if(!new_element) cout << "Error: VoltageSource " << match[1] << " already exists in the circuit" << endl;
                    else if(!exist_control) cout  << "Error: Dependent source has an undefined control element." << endl;
                    else input_type="add_VoltageSource->current";
                }
            }
            else if(regex_search(get,match,pattern[12])){
                if(match[1].str()[0]!='F') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else if(match[3].str()[0]!='N'||match[3].length()==1)cout << "Error: Node " << match[3] << " not found in library" << endl;
                else {
                    bool new_element=1,exist_control=0;
                    circuit.containsElementWithName(new_element,match[1]);
                    for(int i=0;i<circuit.element.size();i++){
                        if(circuit.element[i]->getType()=="V"&&circuit.element[i]->getName()==match[4]){
                            exist_control=1;
                            break;
                        }
                    }
                    if(!new_element) cout << "Error: CurrentSource " << match[1] << " already exists in the circuit" << endl;
                    else if(!exist_control) cout  << "Error: Dependent source has an undefined control element." << endl;
                    else input_type="add_CurrentSource->current";
                }
            }
            else cout << "Error: Syntax error" << endl;
        }
        else if(regex_match(get,regex (R"(\s*delete\b.*)"))){
            regex pattern[]={regex(R"(^\s*delete\s+([RCLDVIEGHF]\w+)\s*$)"),
                             regex(R"(^\s*delete\s+GND\s+(\w+)\s*$)")};
            if(regex_search(get,match,pattern[0])){
                bool exist_element=0;
                circuit.containsElementWithName(exist_element,match[1]);
                if(!exist_element){
                    if(match[1].str()[0]=='R')cout << "Error: Cannot delete resistor; component not found" << endl;
                    else if(match[1].str()[0]=='C')cout << "Error: Cannot delete capacitor; component not found" << endl;
                    else if(match[1].str()[0]=='L')cout << "Error: Cannot delete inductor; component not found" << endl;
                    else if(match[1].str()[0]=='D')cout << "Error: Cannot delete diode; component not found" << endl;
                    else if(match[1].str()[0]=='V'||match[1].str()[0]=='E'||match[1].str()[0]=='H')cout << "Error: Cannot delete VoltageSource; component not found" << endl;
                    else if(match[1].str()[0]=='I'||match[1].str()[0]=='G'||match[1].str()[0]=='F')cout << "Error: Cannot delete CurrentSource; component not found" << endl;
                }
                else input_type="delete_Element";
            }
            else if(regex_search(get,match,pattern[1])){
                if(match[1].str()[0]!='N'||match[1].length()==1)cout << "Error: Node " << match[1] << " not found in library" << endl;
                else{
                    bool exist_element=0;
                    int index=-1;
                    for(int i=0;i<circuit.node.size();i++){
                        if(circuit.node[i].name==match[1]){
                            exist_element=1;
                            index=i;
                            break;
                        }
                    }
                    if(!exist_element) cout << "ERROR: Node " << match[1] << " does not exist in the circuit" << endl;
                    else if(!circuit.node[index].is_ground)cout << "Error: Cannot delete GND; node " << match[1] << " is not GND" << endl;
                    else input_type="delete_GND";
                }
            }
            else cout << "Error: Syntax error" << endl;
        }
        else if(regex_match(get,regex (R"(\s*nodes\s*)"))){
            input_type="nodes";
        }
        else if(regex_match(get,regex (R"(\s*list\b.*)"))){
            regex pattern[]={regex(R"(^\s*list\s*$)"),regex(R"(^\s*list\s*\[\s*(R|C|L|D|V|I|E|G|H|F)\s*\]\s*$)")};
            if(regex_search(get,match,pattern[0])){
                if(circuit.element.size()==0) cout << "No elements have been added" << endl;
                else input_type="total_list";
            }
            else if(regex_search(get,match,pattern[1])){
                if(circuit.element.size()==0) cout << "No elements have been added" << endl;
                else input_type="spacial_list";
            }
            else cout << "Error: Syntax error" << endl;
        }
        else if(regex_match(get,regex (R"(\s*rename\b.*)"))){
            if(regex_search(get,match,regex (R"(^\s*rename\s+node\s+(\w+)\s+(\w+)\s*$)"))){
                if(match[1].str()[0]!='N'||match[1].length()==1)cout << "Error: Node " << match[1] << " not found in library" << endl;
                else if(match[2].str()[0]!='N'||match[2].length()==1)cout << "Error: Node " << match[2] << " not found in library" << endl;
                else{
                    bool exist_node=0,new_node=1;
                    for(int i=0;i<circuit.node.size();i++){
                        if(circuit.node[i].name==match[1])exist_node=1;
                        if(circuit.node[i].name==match[2])new_node=0;
                    }
                    if(!exist_node) cout << "ERROR: Node " << match[1] << " does not exist in the circuit" << endl;
                    else if(!new_node) cout << "ERROR: Node name " << match[2] << " already exists" << endl;
                    else input_type="rename_node";
                }
            }
            else cout << "ERROR: Invalid syntax - correct format:" << endl << "rename node <old_name> <new_name>" << endl;
        }
        else if(regex_match(get,regex (R"(\s*print\b.*)"))){
            regex pattern[]={regex(R"(^\s*print\s+TRAN\s+(?:([-+]?\d+(?:\.\d+)?)\s*s?\s+)?([-+]?\d+(?:\.\d+)?)\s*s?\s+([-+]?\d+(?:\.\d+)?)\s*s?\s+((?:[IV]\(\w+\)\s+)*(?:[IV]\(\w+\)))\s*$)"),
                             regex(R"(^\s*print\s+DC\s+(\w+)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+((?:[IV]\(\w+\)\s+)*(?:[IV]\(\w+\)))\s*$)")};
            int no_error=0;
            if(regex_search(get,match,pattern[0])){
                string s=match[1];
                if(match[1].length()==0)s+="0";
                if(stod(s)<0) cout << "Error: Start time cannot be negative" << endl;
                if(stod(match[2])<=0) cout << "Error: Stop time cannot be zero or negative." << endl;
                else if(stod(match[3])<=0) cout << "Error: Step time cannot be zero or negative." << endl;
                else if(stod(match[3])> stod(match[2])- stod(s)) cout << "Error: Time step is greater than the given range" << endl;
                else no_error=1;
            }
            else if(regex_search(get,match,pattern[1])){
                bool  exist_element=0;
                for(int i=0;i<circuit.element.size();i++){
                    string s =circuit.element[i]->getType();
                    if((s=="V"||s=="I")&&circuit.element[i]->getName()==match[1]){
                        exist_element=1;
                        break;
                    }
                }
                if(!exist_element) cout << " Component " << match[1] << " not found in circuit" << endl;
                else if(stod(match[4])<=0) cout << "Error: Increment cannot be zero or negative." << endl;
                else if(stod(match[3])- stod(match[2])< stod(match[4])) cout << "Error: Increment is greater than the given range" << endl;
                else no_error=2;
            }
            else cout << "Error: Syntax error" << endl;
            if(no_error!=0){
                string input;
                if(no_error==1)input=match[4];
                else input=match[5];
                regex re(R"(([IV])\((\w+)\))");
                bool exist=0;
                for (sregex_iterator it(input.begin(),input.end(),re);it!=sregex_iterator();it++){
                    string it_str=it->str();
                    smatch sm;
                    regex_search(it_str,sm,re);
                    exist=0;
                    if(sm[1]=="I"){
                        circuit.containsElementWithName(exist,sm[2]);
                        if(!exist){
                            cout << "Component " << sm[2] << " not found in circuit" << endl;
                            break;
                        }
                    }
                    else {
                        for(int i=0;i<circuit.node.size();i++){
                            if(circuit.node[i].name==sm[2]){
                                exist=1;
                                break;
                            }
                        }
                        if(!exist){
                            cout << " Node " << sm[2] << " not found in circuit" << endl;
                            break;
                        }
                    }
                }
                if(exist){
                    bool exist_ground=0;
                    for(int i=0;i<circuit.node.size();i++){
                        if(circuit.node[i].is_ground){
                            exist_ground=1;
                            break;
                        }
                    }
                    if(!exist_ground) cout << "Error: No ground node detected in the circuit" << endl;
                    else if(!circuit.isCircuitComplete(0)) cout << "Error: The circuit is discontinuous or has multiple sections." << endl;
                    else {
                        bool exist_control=0;
                        for(int i=0;i<circuit.element.size();i++){
                            if(V_v* p = dynamic_cast<V_v*>(circuit.element[i].get())){
                                exist_control=0;
                                bool helper_bool=0;
                                for(int i=0;i<circuit.node.size();i++){
                                    if(circuit.node[i].name==p->node1.name||circuit.node[i].name==p->node2.name){
                                        if(!helper_bool)helper_bool=1;
                                        else {
                                            exist_control=1;
                                            break;
                                        }
                                    }
                                }
                                if(!exist_control){
                                    cout << "Error: Dependent source " << circuit.element[i]->getName() << " has an undefined control element" << endl;
                                    break;
                                }
                            }
                            else if(I_v* p = dynamic_cast<I_v*>(circuit.element[i].get())){
                                exist_control=0;
                                bool helper_bool=0;
                                for(int i=0;i<circuit.node.size();i++){
                                    if(circuit.node[i].name==p->node1.name||circuit.node[i].name==p->node2.name){
                                        if(!helper_bool)helper_bool=1;
                                        else {
                                            exist_control=1;
                                            break;
                                        }
                                    }
                                }
                                if(!exist_control){
                                    cout << "Error: Dependent source " << circuit.element[i]->getName() << " has an undefined control element" << endl;
                                    break;
                                }
                            }
                            else if(V_i* p = dynamic_cast<V_i*>(circuit.element[i].get())){
                                exist_control=0;
                                circuit.containsElementWithName(exist_control,p->cntr_element);
                                if(!exist_control){
                                    cout << "Error: Dependent source " << circuit.element[i]->getName() << " has an undefined control element" << endl;
                                    break;
                                }
                            }
                            else if(I_i* p = dynamic_cast<I_i*>(circuit.element[i].get())){
                                exist_control=0;
                                circuit.containsElementWithName(exist_control,p->cntr_element);
                                if(!exist_control){
                                    cout << "Error: Dependent source " << circuit.element[i]->getName() << " has an undefined control element" << endl;
                                    break;
                                }
                            }
                        }
                        if(exist_control){
                            if(no_error==1)input_type="print_TRAN";
                            else input_type="print_DC";
                        }
                    }
                }
            }
        }
        else cout << "Error: Syntax error" << endl;
    }
    void handle_input(string get){
        input_type="";
        Error_handling();
        if(input_type.find("add")==0){
            circuit.Add(match,input_type);
            cout << "Element added successfully." << endl;
        }
        else if(input_type.find("delete")==0){
            circuit.Delete(match,input_type);
            cout << "Element successfully deleted." << endl;
        }
        else if(input_type=="nodes"){
            cout << "Available nodes:" << endl;
            for(int i=0;i<circuit.node.size();i++){
                cout <<  circuit.node[i].name ;
                if(i!=circuit.node.size()-1) cout << " , ";
                else cout << endl;
            }
        }
        else if(input_type=="total_list"||input_type=="spacial_list"){
            cout << "Available elements:" << endl;
            for(int i=0;i<circuit.element.size();i++){
                string requested_type,type;
                type=circuit.element[i]->getType();
                if(input_type=="total_list") requested_type=circuit.element[i]->getType();
                else requested_type=match[1];
                if(requested_type==type){
                    if(type=="R"||type=="C"||type=="L"||type=="V"||type=="I"){
                        cout << "Element name : " << circuit.element[i]->getName();
                        cout << " _ Nodes : " << circuit.element[i]->node1.name << ", " << circuit.element[i]->node2.name;
                        cout << " _ Value : " << circuit.element[i]->getValue() << endl;
                    }
                    else if(type=="D"){
                        cout << "Element name : " << circuit.element[i]->getName();
                        cout << " _ Nodes : " << circuit.element[i]->node1.name << ", " << circuit.element[i]->node2.name;
                        cout << " _ Type : ";
                        if(circuit.element[i]->getValue()==0) cout << "D" << endl;
                        else cout << "Z" << endl;
                    }
                    else if(type=="E"||type=="G"){
                        if (V_v* p = dynamic_cast<V_v*>(circuit.element[i].get())){
                            cout << "Element name : " << p->getName();
                            cout << " _ Nodes : " << p->node1.name << ", " << p->node2.name;
                            cout << " _ ControlNodes : " << p->cntr_node1.name << ", " << p->cntr_node2.name;
                            cout << " _ Gain : " << p->getValue() << endl;
                        }
                        if (I_v* p = dynamic_cast<I_v*>(circuit.element[i].get())){
                            cout << "Element name : " << p->getName();
                            cout << " _ Nodes : " << p->node1.name << ", " << p->node2.name;
                            cout << " _ ControlNodes : " << p->cntr_node1.name << ", " << p->cntr_node2.name;
                            cout << " _ Gain : " << p->getValue() << endl;
                        }
                    }
                    else if(type=="H"||type=="F"){
                        if (V_i* p = dynamic_cast<V_i*>(circuit.element[i].get())){
                            cout << "Element name : " << p->getName();
                            cout << " _ Nodes : " << p->node1.name << ", " << p->node2.name;
                            cout << " _ ControlElement : " << p->cntr_element ;
                            cout << " _ Gain : " << p->getValue() << endl;
                        }
                        if (I_i* p = dynamic_cast<I_i*>(circuit.element[i].get())){
                            cout << "Element name : " << p->getName();
                            cout << " _ Nodes : " << p->node1.name << ", " << p->node2.name;
                            cout << " _ ControlElement : " << p->cntr_element ;
                            cout << " _ Gain : " << p->getValue() << endl;
                        }
                    }
                }
                else if((requested_type=="V"&&type=="Vsin")||(requested_type=="I"&&type=="Isin")){
                    if (Vsin* p = dynamic_cast<Vsin*>(circuit.element[i].get())) p->View_profile();
                    if (Isin* p = dynamic_cast<Isin*>(circuit.element[i].get())) p->View_profile();
                }
            }
        }
        else if(input_type=="rename_node"){
            circuit.Rename(match);
            smatch match;
            regex_search(get,match,regex (R"(^\s*rename\s+node\s+(\w+)\s+(\w+)\s*$)"));
            cout << "SUCCESS: Node renamed from " << match[1] << " to " << match[2] << endl;
        }
        else if(input_type=="print_TRAN"){
            vector<vector<string>> wanted_elements={};
            double tstep,tstart,tstop,tspent;
            string input=match[4];
            regex re(R"(([IV])\((\w+)\))");
            for (sregex_iterator it(input.begin(),input.end(),re);it!=sregex_iterator();it++){
                string it_str=it->str();
                smatch sm;
                regex_search(it_str,sm,re);
                wanted_elements.push_back({sm[1],sm[2]});
                if(sm[1]=="I"){
                    for(int i=0;i<circuit.element.size();i++){
                        if(circuit.element[i]->getName()==sm[2]){
                            wanted_elements[wanted_elements.size()-1].push_back(to_string(i));
                        }
                    }
                }
                else {
                    for(int i=0;i<circuit.node.size();i++){
                        if(circuit.node[i].name==sm[2]){
                            wanted_elements[wanted_elements.size()-1].push_back(to_string(i));
                        }
                    }
                }
            }
            tstep= stod(match[3]);
            tstop= stod(match[2]);
            string s_tstart=match[1];
            if(s_tstart.length()==0)s_tstart.push_back('0');
            tstart= stod(s_tstart);
            int number=(tstop-tstart)/tstep;
            tspent=tstart;
            for(int i=0;i<=number;i++){
                tspent+=double(i)*tstep;
                if(circuit.Print_TRAN(match,tspent)){
                    cout << tspent <<  " seconds have passed since the start of the circuit :" << endl;
                    for(int j=0;j<wanted_elements.size();j++){
                        cout << wanted_elements[j][1] << " : " << wanted_elements[j][0] << " = ";
                        if(wanted_elements[j][0]=="V"){
                            if(tspent==0.0)cout << circuit.matrixSolve.Primary_TRAN.Answer[stoi(wanted_elements[j][2])]<< " volt." << endl;
                            else cout << circuit.matrixSolve.TRAN.Answer[stoi(wanted_elements[j][2])]<< " volt." << endl;
                        }
                        else {
                            int index= stoi(wanted_elements[j][2]);
                            if(tspent==0.0)cout << fixed << setprecision(4) <<circuit.element_current_shower(index,circuit.matrixSolve.Primary_TRAN.Answer,tspent,tstep,"TRAN");
                            else cout << fixed << setprecision(4) <<circuit.element_current_shower(index,circuit.matrixSolve.TRAN.Answer,tspent,tstep,"TRAN");
                            cout << " amp." << endl;
                        }
                    }
                }
                else {
                    cout << "Error: The electrical circuit cannot be solved." << endl;
                    break;
                }
            }
        }
        else if(input_type=="print_DC"){
            vector<vector<string>> wanted_elements={};
            double StartValue,EndValue,Increment;
            for(int i=0;i<circuit.element.size();i++){
                if(circuit.element[i]->getName()==match[1]){
                    if(circuit.element[i]->getType()=="V")wanted_elements.push_back({"V",match[1], to_string(i)});
                    else wanted_elements.push_back({"I",match[1], to_string(i)});
                }
            }
            string input=match[5];
            regex re(R"(([IV])\((\w+)\))");
            for (sregex_iterator it(input.begin(),input.end(),re);it!=sregex_iterator();it++){
                string it_str=it->str();
                smatch sm;
                regex_search(it_str,sm,re);
                wanted_elements.push_back({sm[1],sm[2]});
                if(sm[1]=="I"){
                    for(int i=0;i<circuit.element.size();i++){
                        if(circuit.element[i]->getName()==sm[2]){
                            wanted_elements[wanted_elements.size()-1].push_back(to_string(i));
                        }
                    }
                }
                else {
                    for(int i=0;i<circuit.node.size();i++){
                        if(circuit.node[i].name==sm[2]){
                            wanted_elements[wanted_elements.size()-1].push_back(to_string(i));
                        }
                    }
                }
            }
            StartValue= stod(match[2]);
            EndValue= stod(match[3]);
            Increment= stod(match[4]);
            int number=(EndValue-StartValue)/Increment;
            for(int i=0;i<=number;i++){
                StartValue+=double(i)*Increment;
                if(circuit.Print_DC(StartValue, stoi(wanted_elements[0][2]))){
                    if(wanted_elements[0][0]=="V") cout << "VoltageSource : "<< wanted_elements[0][1] << " _ Voltage : " << StartValue << " volt. DC :"<< endl;
                    else cout << "CurrentSource : " << wanted_elements[0][1] << " _ current : " << StartValue << " amp. DC :"<< endl;;
                    for(int j=1;j<wanted_elements.size();j++){
                        cout << wanted_elements[j][1] << " : " << wanted_elements[j][0] << " = ";
                        if(wanted_elements[j][0]=="V"){
                            cout << circuit.matrixSolve.DC.Answer[stoi(wanted_elements[j][2])]<< " volt." << endl;
                        }
                        else {
                            int index= stoi(wanted_elements[j][2]);
                            circuit.element_current_shower(index,circuit.matrixSolve.DC.Answer,0,0,"DC");
                            cout << " amp." << endl;
                        }
                    }
                }
                else {
                    cout << "Error: The electrical circuit cannot be solved." << endl;
                    break;
                }
            }
        }
    }
public:
    void Run(){
        while (true){
            getline(cin,get);
            if(get=="Exit"){
                cout << "Good luck!" << endl;
                break;
            }
            handle_input(get);
        }
    }
};
int main() {
    View view;
    view.Run();
    return 0;
}