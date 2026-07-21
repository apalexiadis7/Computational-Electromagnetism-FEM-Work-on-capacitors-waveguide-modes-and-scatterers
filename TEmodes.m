function[]=TEmodes(refinements)
% ορισμός μιας απλής κυκλικής γεωμετρίας
gd=[1, 0, 0, 1e-2]';
d1=decsg(gd);
[p,e,t] = initmesh(d1);
% σε περίπτωση που θέλουμε επιπλέον πύκνωση του πλέγματος
if refinements~=0
for i=1:refinements
[p,e,t]=refinemesh(d1,p,e,t);
end
end
Nn=size(p,2);
Nd=size(e,2);
Ne=size(t,2);
% αρχικοποιούμε με μηδενικά το διάνυσμα των Ηz σε όλον τον χώρο. Έχουμε
% οριακές συνθήκες Neumann οπότε όλοι οι κόμβοι είναι άγνωστοι
H0 = zeros(Nn,6);
figure;
pdeplot(p,e,t);
%αρχικοποίηση των ολικών αραιών πινάκων ακαμψίας και μάζης
S = spalloc(Nn,Nn,7*Nn);
T = spalloc(Nn,Nn,7*Nn);
% Ακολουθεί ο πυρήνας του κώδικα, όπως αναγράφεται στις διαφάνειες
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
%τοπικός πίνακας ακαμψίας
Te=zeros(3,3);
%τοπικός πίνακας μάζας
for i=1:3
for j=1:3
Se(i,j) = (b(i)*b(j)+c(i)*c(j))*Ae;
if i==j
Te(i,j) = Ae/6;
else
Te(i,j) = Ae/12;
end
S((n(i)),(n(j))) = S((n(i)),(n(j))) + Se(i,j);
T((n(i)),(n(j))) = T((n(i)),(n(j))) + Te(i,j);
end
end
end
% επίλυση του προβλήματος ιδιοτιμών, λαμβάνοντας τις 10 μικρότερες
sigma='smallestabs';
% με αυτή την έκφραση της συνάρτησης, μας επιστρέφεται ένα διάνυσμα που
% περιέχει τις ιδιοτιμές τις οποίες ακριβώς μετά πλοτάρουμε
d=eigs(S,T,10,sigma);
figure;
plot(d,'o');
% αυτή η έκφραση της συνάρτησης επιστρέφει έναν τετραγωνικό πίνακα D με τις
% ιδιοτιμές του στην κύρια διαγώνιο και έναν πίνακα V με στήλες τα
% αντίστοιχα ιδιοδιανύσματα
[V,D]=eigs(S,T,10,sigma);
% εδώ γίνεται επιλογή των καταλλήλων ιδιοτιμών/ιδιοδιανυσμάτων
V1=zeros(Nn,6);
for i=1:Nn
V1(i,1)=V(i,1);
V1(i,2)=V(i,2);
V1(i,3)=V(i,4);
V1(i,4)=V(i,6);
V1(i,5)=V(i,7);
V1(i,6)=V(i,9);
end
% αποτυπώνονται οι ρυθμοί από τα ιδιοδιανύσματα
for i=1:6
for j=1:Nn
H0(j,i)=V1(j,i);
end
end
% εδώ παράγονται έξι διαφορετικά διανύσματα, ουσιαστικά οι στήλες του H0,
% όπου το καθένα αντιπροσωπέυει έναν συγκεκριμένο ρυθμό
X1=zeros(Nn,1);
X2=zeros(Nn,1);
X3=zeros(Nn,1);
X4=zeros(Nn,1);
X5=zeros(Nn,1);
X6=zeros(Nn,1);
for i=1:Nn
X1(i,1)=H0(i,1);
X2(i,1)=H0(i,2);
X3(i,1)=H0(i,3);
X4(i,1)=H0(i,4);
X5(i,1)=H0(i,5);
X6(i,1)=H0(i,6);
end
% απεικόνιση της z συνιστώσας
for i=1:6
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',H0(1:Nn,i),'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
end
% βρίσκουμε το gradient, εφόσον το εγκάρσιο πεδίο δίνεται από αυτό. Συγκεκριμένα η
% μερική παράγωγος κατά y είναι ανάλογη του Εx ενώ η μερική παράγωγος κατά
% x έιναι ανάλογη του Ey
[gradx1,grady1] = pdegrad(p,t,X1);
[gradx2,grady2] = pdegrad(p,t,X2);
[gradx3,grady3] = pdegrad(p,t,X3);
[gradx4,grady4] = pdegrad(p,t,X4);
[gradx5,grady5] = pdegrad(p,t,X5);
[gradx6,grady6] = pdegrad(p,t,X6);
% μας αφορά το μέτρο του πεδίου, το οποίο προκύπτει από τη ρίζα των
% τετραγώνων των Ex & Ey
u1=zeros(Ne,1);
u2=zeros(Ne,1);
u3=zeros(Ne,1);
u4=zeros(Ne,1);
u5=zeros(Ne,1);
u6=zeros(Ne,1);
for i=1:Ne
u1(i,1)=sqrt(((gradx1(i))^2)+((grady1(i))^2));
u2(i,1)=sqrt(((gradx2(i))^2)+((grady2(i))^2));
u3(i,1)=sqrt(((gradx3(i))^2)+((grady3(i))^2));
u4(i,1)=sqrt(((gradx4(i))^2)+((grady4(i))^2));
u5(i,1)=sqrt(((gradx5(i))^2)+((grady5(i))^2));
u6(i,1)=sqrt(((gradx6(i))^2)+((grady6(i))^2));
end
% αναπαράσταση των επιθυμητών ρυθμών
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',gradx1,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u2,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u3,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u4,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u5,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u6,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
end
