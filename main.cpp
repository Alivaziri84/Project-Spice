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
class Node{
public:
    string name;
    int index;
    bool is_ground=0;
};
class Element{
protected:
    string type;
    string name;
    double value;
public:
    Node node1,node2;
    Element(string t,string n, double v) : type(t),name(n),value(v){}
    virtual void Add(vector<vector<double>> & G,vector<double> I)=0;
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
    Resistor(string t,string n, double v) : Element(t,n,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{}
};
class Capacitor : public Element{
public:
    Capacitor(string t,string n, double v) : Element(t,n,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{}
};
class Inductor : public Element{
public:
    Inductor(string t,string n, double v) : Element(t,n,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{
    }
};
class VoltageSource : public Element{
public:
    VoltageSource(string t,string n, double v) : Element(t,n,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{}
};
class CurrentSource : public Element{
public:
    CurrentSource(string t,string n, double v) : Element(t,n,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{}
};
class Diode : public Element{
public:
    Diode(string t,string n, double v) : Element(t,n,v){}
    void Add(vector<vector<double>> & left,vector<double> right) override{}
};
class Circuit{
private:
    Matrix_solve matrixSolve;
    /*
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
     */
public:
    void containsElementWithName(bool &b, string s){
        for(int i=0;i<element.size();i++){
            if(element[i]->getName()==s){
                b=!b;
                break;
            }
        }
    }
    vector<unique_ptr<Element>> element={};
    vector<Node> node={};
    void Add(smatch match,string input_type){
        if(input_type=="add_Resistor"){}
        else if(input_type=="add_Capacitor"){}
        else if(input_type=="add_Inductor"){}
        else if(input_type=="add_Diode"){}
        else if(input_type=="add_GND"){}
        else if(input_type=="add_VoltageSource"){}
        else if(input_type=="add_CurrentSource"){}
        else if(input_type=="add_VoltageSource_sin"){}
        else if(input_type=="add_CurrentSource_sin"){}
        else if(input_type=="add_VoltageSource->voltage"){}
        else if(input_type=="add_CurrentSource->voltage"){}
        else if(input_type=="add_VoltageSource->current"){}
        else if(input_type=="add_CurrentSource->current"){}
    }
    void Delete(smatch match,string input_type){
        if(input_type=="delete_Element"){}
        else if(input_type=="delete_GND"){}
    }
    void Rename(smatch match){

    }
    void Print(smatch match,string input_type){
        if(input_type=="TRAN_print"){}
        else if(input_type=="DC_print"){}
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
                             regex(R"(^\s*add\s+([^I]\w+)\s+(\w+)\s+(\w+)\s+SIN\(\s*([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*\)\s*$)"),
                             regex(R"(^\s*add\s+([^V]\w+)\s+(\w+)\s+(\w+)\s+SIN\(\s*([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s+([-+]?\d+(?:\.\d+)?(?:Meg|[kMunm]|e[-+]?\d+(?:\.\d+)?)?)\s*\)\s*$)"),
                             regex(R"(^\s*add\s+([^G]\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)"),//E
                             regex(R"(^\s*add\s+([^E]\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)"),//G
                             regex(R"(^\s*add\s+([^F]\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)"),//H
                             regex(R"(^\s*add\s+([^H]\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+([-+]?\d+(?:\.\d+)?)\s*$)")//F
            };
            if(regex_search(get,match,pattern[0])){
                if(match[1].str()[0]!='R') cout << "Error: Element " << match[1] << " not found in library" << endl;
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
                else if(stod(match[4])<=0) cout << "Error: Inductance cannot be zero or negative" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: inductor " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_Inductor";
                }
            }
            else if(regex_search(get,match,pattern[3])){
                if(match[4]!="D"&&match[4]!="Z") cout << "Error: Model " << match[4] << " not found in library" << endl;
                else if(match[1].str()[0]!='D') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: diode " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_Diode";
                }
            }
            else if(regex_search(get,match,pattern[4])){
                if(match[1]!="GND") cout << "Error: Element GND not found in library" << endl;
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
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: VoltageSource " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_VoltageSource";
                }
            }
            else if(regex_search(get,match,pattern[6])){
                if(match[1].str()[0]!='I') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: CurrentSource " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_CurrentSource";
                }
            }
            else if(regex_search(get,match,pattern[7])){
                if(match[1].str()[0]!='V') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: VoltageSource " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_VoltageSource_sin";
                }
            }
            else if(regex_search(get,match,pattern[8])){
                if(match[1].str()[0]!='I') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else {
                    bool new_element=1;
                    circuit.containsElementWithName(new_element,match[1]);
                    if(!new_element) cout << "Error: CurrentSource " << match[1] << " already exists in the circuit" << endl;
                    else input_type="add_CurrentSource_sin";
                }
            }
            else if(regex_search(get,match,pattern[9])){
                if(match[1].str()[0]!='E') cout << "Error: Element " << match[1] << " not found in library" << endl;
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
                else {
                    bool new_element=1,exist_control=0;
                    circuit.containsElementWithName(new_element,match[1]);
                    circuit.containsElementWithName(exist_control,match[4]);
                    if(!new_element) cout << "Error: VoltageSource " << match[1] << " already exists in the circuit" << endl;
                    else if(!exist_control) cout  << "Error: Dependent source has an undefined control element." << endl;
                    else input_type="add_VoltageSource->current";
                }
            }
            else if(regex_search(get,match,pattern[12])){
                if(match[1].str()[0]!='F') cout << "Error: Element " << match[1] << " not found in library" << endl;
                else {
                    bool new_element=1,exist_control=0;
                    circuit.containsElementWithName(new_element,match[1]);
                    circuit.containsElementWithName(exist_control,match[4]);
                    if(!new_element) cout << "Error: VoltageSource " << match[1] << " already exists in the circuit" << endl;
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
                bool exist_node=0,new_node=1;
                for(int i=0;i<circuit.node.size();i++){
                    if(circuit.node[i].name==match[1])exist_node=1;
                    if(circuit.node[i].name==match[2])new_node=0;
                }
                if(!exist_node) cout << "ERROR: Node " << match[1] << " does not exist in the circuit" << endl;
                else if(!new_node) cout << "ERROR: Node name " << match[2] << " already exists" << endl;
                else input_type="rename_node";
            }
            else cout << "ERROR: Invalid syntax - correct format:" << endl << "rename node <old_name> <new_name>" << endl;
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
    void handle_input(string get){
        input_type="";
        Error_handling();
        if(input_type.find("add")==0){
            circuit.Add(match,input_type);
        }
        else if(input_type.find("delete")==0){
            circuit.Delete(match,input_type);
        }
        else if(input_type=="nodes"){
            cout << "Available nodes:" << endl;
            for(int i=0;i<circuit.node.size();i++){
                cout <<  circuit.node[i].name << endl;
                if(i!=circuit.node.size()-1) cout << " , ";
                else cout << endl;
            }
        }
        else if(input_type.find("list")!=-1){
            if(input_type=="total_list"){
                cout << "Available elements:" << endl;
                for(int i=0;i<circuit.element.size();i++){
                    string type=circuit.element[i]->getType();
                    //RCLDVIEGHF
                    if(type=="R"||type=="C"||type=="L"||type=="V"||type=="I"){
                        cout << "Element name : " << circuit.element[i]->getName();
                        cout << " _ Nodes : " << circuit.element[i]->node1.name << ", " << circuit.element[i]->node2.name;
                        cout << " _ Value : " << circuit.element[i]->getValue() << endl;
                    }
                }
            }
            else if(input_type=="spacial_list"){}
        }
        else if(input_type=="rename_node"){
            circuit.Rename(match);
            smatch match;
            regex_search(get,match,regex (R"(^\s*rename\s+node\s+(\w+)\s+(\w+)\s*$)"));
            cout << "SUCCESS: Node renamed from " << match[1] << " to " << match[2] << endl;
        }
        else if(input_type.find("print")==0){
            //circuit.Print(get);
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