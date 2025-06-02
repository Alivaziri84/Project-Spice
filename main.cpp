#include <iostream>
#include <regex>
#include <vector>
#include <cmath>
using namespace std;
class Matrix_solve{
public:
    vector<vector<double>> left_matrix={};
    vector<double> answer={};
    vector<double> right_matrix={};
    vector<vector<double>> l_left={};
    vector<vector<double>> u_left={};
    void LU_setter(){
        int size=left_matrix.size();
        for(int i=0;i<size;i++){
            for(int j=i;j<size;j++){
                if(fabs(left_matrix[i][i])<fabs(left_matrix[j][i]))swap(left_matrix[i],left_matrix[j]);
            }
        }
        l_left={};
        u_left={};
        for(int i=0;i<size;i++){
            l_left.push_back({});
            u_left.push_back({});
            for(int j=0;j<size;j++){
                if(j==i) l_left[i].push_back(1.0);
                else l_left[i].push_back(0.0);
                u_left[i].push_back(0.0);
            }
        }
        for(int i=0;i<size;i++){
            for(int j=i;j<size;j++){
                double sumu=0.0,suml=0.0;
                for(int k=0;k<size;k++){
                    sumu+=l_left[i][k]*u_left[k][j];
                }
                u_left[i][j]=left_matrix[i][j]-sumu;
                if(j!=i){
                    for(int k=0;k<size;k++){
                        suml+=l_left[j][k]*u_left[k][i];
                    }
                    l_left[j][i]=(left_matrix[j][i]-suml)/u_left[i][i];
                }
            }
        }
    }
    void Solve(){
        if(0){
            vector<vector<double>> v={};
            for(int i=0;i<left_matrix.size();i++){
                v.push_back(left_matrix[i]);
                v[i].push_back(right_matrix[i]);
            }
            int size=v.size();
            for(int i=0;i<size-1;i++){
                if(v[i][i]<1e-9){
                    for(int j=i;j<size;j++){
                        if(v[j][i]>1e-9){
                            swap(v[i],v[j]);
                            break;
                        }
                    }
                }
                for(int j=i+1;j<size;j++){
                    for(int k=0;k<size+1;k++){
                        if(k!=i) v[j][k]-=v[i][k]*v[j][i]/v[i][i];
                    }
                    v[j][i]=0.0;
                }
            }
            for(int i=size-2;i>=0;i--){
                for(int j=i+1;j<size;j++){
                    v[i][size]-=v[j][size]*v[i][j]/v[j][j];
                    v[i][j]=0.0;
                }
            }
            answer={};
            for(int i=0;i<size;i++) answer.push_back(v[i][size]/v[i][i]);
        }
        else{
            int size=right_matrix.size();
            vector<vector<double>> v={};
            for(int i=0;i<size;i++){
                v.push_back(l_left[i]);
                v[i].push_back(right_matrix[i]);
            }
            for(int i=1;i<size;i++){
                for(int j=0;j<i;j++){
                    v[i][size]-=v[i][j]*v[j][size];
                }
            }
            for(int i=0;i<size;i++){
                for(int j=0;j<size;j++){
                    v[i][j]=u_left[i][j];
                }
            }
            for(int i=size-2;i>=0;i--){
                for(int j=i+1;j<size;j++){
                    v[i][size]-=v[j][size]*v[i][j]/v[j][j];
                }
            }
            answer={};
            for(int i=0;i<size;i++) answer.push_back(v[i][size]/v[i][i]);
        }
    }
};
class Element{
protected:
    string type;
    string name;
    int node1,node2;
    double value;
public:
    Element(string t,string n,int n1,int n2, double v) : type(t),name(n),node1(n1),node2(n2),value(v){}
    virtual void Add(vector<vector<double>> & G,vector<double> I)=0;
    int getIndex1(){
        return node1;
    };
    int getIndex2(){
        return node2;
    };
    string getName(){
        return name;
    }
    void setValue(double v){
        if(v>0)value=v;
        else{}
    }
    double getValue(){
        return value;
    };
    string getType(){
        return type;
    }
};
class Resistor : public Element{
public:
    Resistor(string t,string n,int n1,int n2, double v) : Element(t,n,n1,n2,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{}
};
class Capacitor : public Element{
public:
    Capacitor(string t,string n,int n1,int n2, double v) : Element(t,n,n1,n2,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{}
};
class Inductor : public Element{
public:
    Inductor(string t,string n,int n1,int n2, double v) : Element(t,n,n1,n2,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{
    }
};
class VoltageSource : public Element{
public:
    VoltageSource(string t,string n,int n1,int n2, double v) : Element(t,n,n1,n2,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{}
};
class CurrentSource : public Element{
public:
    CurrentSource(string t,string n,int n1,int n2, double v) : Element(t,n,n1,n2,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{}
};
class Diode : public Element{
public:
    Diode(string t,string n,int n1,int n2, double v) : Element(t,n,n1,n2,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{}
};
class Node{
public:
    string name;
    int index;
    vector<int> NeighbourNodes={};
};
class Circuit{
private:
    vector<Element*> element={};
    vector<Resistor> resistor={};
    vector<Capacitor> capacitor={};
    vector<Inductor> inductor={};
    vector<VoltageSource> voltageSource={};
    vector<CurrentSource> currentSource={};
    vector<Diode> diode={};
    vector<Node> node={};
    Matrix_solve matrixSolve;
    bool isCircuitComplete(bool first,int index){
        static vector<bool> v;
        if(first){
            for(int i=0;i<node.size();i++){
                if(node[i].NeighbourNodes.size()==1)return false;
            }
            v.clear();
            for(int i=0;i<node.size();i++){
                v.push_back(false);
            }
        }
        v[index]=true;
        for(int i=0;i<node[index].NeighbourNodes.size();i++){
            int n=node[index].NeighbourNodes[i];
            if (!v[n]) isCircuitComplete(false,n);
        }
        if(first){
            for(int i=0;i<v.size();i++){
                if(!v[i]) return false;
            }
            return true;
        }
        return false;
    }
public:
    string input_type="null";
    void Erorr_handling(string get){
        smatch match;
        if(regex_match(get,regex (R"(^\s*add\b.*)"))) {
            regex pattern[]={regex(R"(^\s*add\s+(R\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*$)"),
                             regex(R"(^\s*add\s+(C\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*$)"),
                             regex(R"(^\s*add\s+(L\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*$)"),
                             regex(R"(^\s*add\s+(D\w+)\s+(\w+)\s+(\w+)\s+([DZ])\s*$)"),
                             regex(R"(^\s*add\s+GND\s+(\w+)\s*$)"),
                             regex(R"(^\s*add\s+VoltageSource(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*$)"),
                             regex(R"(^\s*add\s+CurrentSource(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*$)"),
                             regex(R"(^\s*add\s+V(\w+)\s+(\w+)\s+(\w+)\s+SIN\(\s*([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*\)\s*$)"),
                             regex(R"(^\s*add\s+E(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)"),
                             regex(R"(^\s*add\s+G(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)"),
                             regex(R"(^\s*add\s+H(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)"),
                             regex(R"(^\s*add\s+F(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)")
            };
            if(regex_search(get,match,pattern[0])){}
            else if(regex_search(get,match,pattern[1])){}
            else if(regex_search(get,match,pattern[2])){}
            else if(regex_search(get,match,pattern[3])){}
            else if(regex_search(get,match,pattern[4])){}
            else if(regex_search(get,match,pattern[5])){}
            else if(regex_search(get,match,pattern[6])){}
            else if(regex_search(get,match,pattern[7])){}
            else if(regex_search(get,match,pattern[8])){}
            else if(regex_search(get,match,pattern[9])){}
            else if(regex_search(get,match,pattern[10])){}
            else if(regex_search(get,match,pattern[11])){}
            else cout << "Error: Syntax error" << endl;
        }
        else if(regex_match(get,regex (R"(\s*delete\b.*)"))){
            regex pattern[]={regex(R"(^\s*delete\s+(R\w+)\s*$)"),
                             regex(R"(^\s*delete\s+(C\w+)\s*$)"),
                             regex(R"(^\s*delete\s+(L\w+)\s*$)"),
                             regex(R"(^\s*delete\s+(D\w+)\s*$)"),
                             regex(R"(^\s*delete\s+GND\s+(\w+)\s*$)"),
                             regex(R"(^\s*delete\s+VoltageSource(\w+)\s*$)"),
                             regex(R"(^\s*delete\s+CurrentSource(\w+)\s*$)")};
            if(regex_search(get,match,pattern[0])){}
            else if(regex_search(get,match,pattern[1])){}
            else if(regex_search(get,match,pattern[2])){}
            else if(regex_search(get,match,pattern[3])){}
            else if(regex_search(get,match,pattern[4])){}
            else if(regex_search(get,match,pattern[5])){}
            else if(regex_search(get,match,pattern[6])){}
            else cout << "Error: Syntax error" << endl;
        }
        else if(regex_match(get,regex (R"(\s*nodes\s*)"))){
            input_type="nodes";
        }
        else if(regex_match(get,regex (R"(\s*list\b.*)"))){
            regex pattern[]={regex(R"(^\s*list\s*$)"),regex(R"(^\s*list\s*\[\s*(R|C|L|D|VoltageSource|CurrentSource)\s*\]\s*$)")};
            if(regex_search(get,match,pattern[0])){}
            else if(regex_search(get,match,pattern[1])){}
            else cout << "Error: Syntax error" << endl;
        }
        else if(regex_match(get,regex (R"(\s*rename\b.*)"))){
            if(regex_search(get,match,regex (R"(^\s*rename\s+node\s+(\w+)\s+(\w+)\s*$)"))){}
            else cout <<"Error: Syntax error" << endl;
        }
        else if(regex_match(get,regex (R"(\s*print\b.*)"))){
            regex pattern[]={regex(R"(^\s*print\s+TRAN$)"),
                             regex(R"()")};
            if(regex_search(get,match,pattern[0])){}
            else if(regex_search(get,match,pattern[1])){}
            else cout << "Error: Syntax error" << endl;
        }
        else cout << "Error: Syntax error" << endl;
    }
    void Add(string get){
        if(input_type=="add_Resistor"){}
        else if(input_type=="add_Capacitor"){}
        else if(input_type=="add_Inductor"){}
        else if(input_type=="add_Diode"){}
        else if(input_type=="add_GND"){}
        else if(input_type=="add_VoltageSource"){}
        else if(input_type=="add_CurrentSource"){}
        else if(input_type=="add_VoltageSource_sin"){}
        else if(input_type=="add_VoltageSource->voltage"){}
        else if(input_type=="add_VoltageSource->current"){}
        else if(input_type=="add_CurrentSource->voltage"){}
        else if(input_type=="add_CurrentSource->current"){}
    }
    void Delete(string get){
        if(input_type=="delete_Resistor"){}
        else if(input_type=="delete_Capacitor"){}
        else if(input_type=="delete_Inductor"){}
        else if(input_type=="delete_Diode"){}
        else if(input_type=="delete_GND"){}
        else if(input_type=="delete_VoltageSource"){}
        else if(input_type=="delete_CurrentSource"){}
    }
    void Nodes(){
        if(input_type=="nodes"){}
    }
    void List(string get){
        if(input_type=="total_list"){}
        else if(input_type=="spacial_list"){}
    }
    void Rename(string get){
        if(input_type=="rename"){}
    }
    void Print(){
        if(input_type=="TRAN_print"){}
        else if(input_type=="DC_print"){}
    }
};

class View{
private:
    Circuit circuit;
public:
    void Run(){
        string get="";
        while (true){
            getline(cin,get);
            if(get=="Exit"){
                cout << "Good luck!" << endl;
                break;
            }
            circuit.Erorr_handling(get);
            if(circuit.input_type.find("add")==0){
                circuit.Add(get);
            }
            else if(circuit.input_type.find("delete")==0){
                circuit.Delete(get);
            }
            else if(circuit.input_type=="nodes"){
                circuit.Nodes();
            }
            else if(circuit.input_type.find("list")==0){
                circuit.List(get);
            }
            else if(circuit.input_type.find("rename")==0){
                circuit.Rename(get);
            }
            else if(circuit.input_type.find("print")==0){
                circuit.Print();
            }
        }
    }
};
int main() {
    View view;
    view.Run();
    return 0;
}