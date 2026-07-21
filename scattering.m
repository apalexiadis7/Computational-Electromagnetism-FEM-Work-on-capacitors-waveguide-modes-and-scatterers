function[]=scattering(a,R,refinements)
% ορισμός δυο κυκλικών γεωμετριών
R1=[1, 0, 0, R]';
R2=[1, 0, 0, a]';
gd=[R1,R2];
ns=(char('R1','R2'))';
d1=decsg(gd,'R1-R2',ns);
[p,e,t] = initmesh(d1);
% σε περίπτωση που θέλουμε επιπλέον πύκνωση του πλέγματος
if refinements~=0
for i=1:refinements
[p,e,t]=refinemesh(d1,p,e,t);
end
end
pdeplot(p,e,t);
e0=8.854e-12;
% διηλεκτρική σταθερά του κενού
m0=4*pi*1e-7;
% μαγνητική διαπερατότητα του κενού
c0=3e8;
% ταχύτητα φωτός
f=300e6;
% συχνότητα
k0=2*pi*f/c0;
% κυματικός αριθμός
Nn=size(p,2);
Nd=size(e,2);
Ne=size(t,2);
node_id = ones(Nn,1);
Ei = zeros(Nn,1);
% το (γνωστό) προσπίπτον πεδίο
Es = zeros(Nn,1);
% το σκεδαζόμενο πεδίο
% υπολογισμός του γνωστού πεδίου
for i=1:Nn
Ei(i)=(exp(1))^(-1i*k0*p(1,i));
end
% Οι κόμβοι που βρίσκονται στο εσωτερικό όριο είναι γνωστοί
for i=1:Nd
if sqrt(((p(1,e(1,i)))^2)+((p(2,e(1,i)))^2))<(a+R)/2
% θεωρούμε τη συνθήκη βάσει του μέσου των δύο ακτινών
Es(e(1,i))=-Ei(e(1,i));
Es(e(2,i))=-Ei(e(2,i));
node_id(e(1,i))=0;
node_id(e(2,i))=0;
end
end
% επαναρίθμηση αγνώστων κόμβων με το index
counter=0;
for id=1:Nn
if node_id(id)==1
counter=counter+1;
index(id)=counter;
end
end
Nf=counter;
% Ακολουθεί ο πυρήνας του κώδικα, όπως αναγράφεται στις διαφάνειες
S = spalloc(Nf,Nf,7*Nf);
T = spalloc(Nf,Nf,7*Nf);
% ο επιπλέον πίνακας του αριστερού μέλους. Θέλουμε ABC 1ης τάξης, επομένως
% μας αφορά μόνον ο πίνακας μάζης πάνω στο εξωτερικό όριο
Tc = spalloc(Nf,Nf,7*Nf);
A = spalloc(Nf,Nf,7*Nf);
% το δεξί μέλος
B1 = zeros(Nf,1);
B2 = zeros(Nf,1);
B = zeros(Nf,1);
for ie = 1:Ne
n(1:3) = t(1:3,ie);
x(1:3) = p(1,n(1:3));
y(1:3) = p(2,n(1:3));
De = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);
Ae = abs(De/2);
b(1) = (y(2)-y(3))/De;
c(1) = (x(3)-x(2))/De;
b(2) = (y(3)-y(1))/De;
c(2) = (x(1)-x(3))/De;
b(3) = (y(1)-y(2))/De;
c(3) = (x(2)-x(1))/De;
Se=zeros(3,3);
Te=zeros(3,3);
for i=1:3
for j=1:3
Se(i,j) = (m0^(-1))*(b(i)*b(j)+c(i)*c(j))*Ae;
if i==j
Te(i,j) = Ae/6;
else
Te(i,j) = Ae/12;
end
if (node_id(n(i))==1)
if (node_id(n(j))==1)
S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j);
T(index(n(i)),index(n(j))) = T(index(n(i)),index(n(j))) + Te(i,j);
else
B1(index(n(i))) = B1(index(n(i))) - Se(i,j)*Es(n(j));
B2(index(n(i))) = B2(index(n(i))) + Te(i,j)*Es(n(j));
end
end
end
end
end
% σχηματισμός Tc που αφορά το εξωτερικό όριο
for id=1:Nd
if sqrt(((p(1,e(1,id)))^2)+((p(2,e(1,id)))^2))>(a+R)/2
n(1:2) = e(1:2,id);
x(1:2) = p(1,n(1:2));
y(1:2) = p(2,n(1:2));
Le = sqrt(((x(1)-x(2))^2)+((y(1)-y(2))^2));
Te = zeros(2,2);
for i=1:2
for j=1:2
if i==j
Te(i,j)=Le/3;
else
Te(i,j)=Le/6;
end
Tc(index(n(i)),index(n(j))) = T(index(n(i)),index(n(j))) + Te(i,j);
end
end
end
end
% υπολογισμός των συντελεστών για το εξωτερικό όριο
a_R=(-1i*k0)-(0.5/R)+(0.125/(1i*k0*(R^2)))+(0.125/((k0^2)*(R^3)));
% υπολογισμός των πινάκων
A = S-(((2*pi*f)^2)*e0*T)-((a_R/m0)*Tc);
B = B1+((2*pi*f)^2)*e0*B2;
X = A\B;
for i=1:Nn
if Es(i)==0
Es(i)=X(index(i));
end
end
% υπολογισμός συνολικού πεδίου
Etot=zeros(Nn,1);
for i=1:Nn
Etot(i)=Es(i)+Ei(i);
end
% απεικόνιση μέτρου και φάσης του ΣΥΝΟΛΙΚΟΥ πεδίου
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',abs(Etot),'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',angle(Etot),'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
end
Ακολουθού
