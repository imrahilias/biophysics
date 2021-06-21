data = dlmread( 'crall_013_blue.csv' )
d2 = data(2:end,:);
d2(:,2) = sqrt(d2(:,2));
d2(:,3) = sqrt(d2(:,3));
d2(:,4) = sqrt(d2(:,4));
dlmwrite( 'crall_013_blue2.csv', d2 )

data = dlmread( 'crall_017_blue.csv' )
d2 = data(2:end,:);
d2(:,2) = sqrt(d2(:,2));
d2(:,3) = sqrt(d2(:,3));
d2(:,4) = sqrt(d2(:,4));
dlmwrite( 'crall_017_blue2.csv', d2 )

data = dlmread( 'crall_019_blue.csv' )
d2 = data(2:end,:);
d2(:,2) = sqrt(d2(:,2));
d2(:,3) = sqrt(d2(:,3));
d2(:,4) = sqrt(d2(:,4));
dlmwrite( 'crall_019_blue2.csv', d2 )

data = dlmread( 'crall_013_red.csv' )
d2 = data(2:end,:);
d2(:,2) = sqrt(d2(:,2));
d2(:,3) = sqrt(d2(:,3));
d2(:,4) = sqrt(d2(:,4));
dlmwrite( 'crall_013_red2.csv', d2 )

data = dlmread( 'crall_017_red.csv' )
d2 = data(2:end,:);
d2(:,2) = sqrt(d2(:,2));
d2(:,3) = sqrt(d2(:,3));
d2(:,4) = sqrt(d2(:,4));
dlmwrite( 'crall_017_red2.csv', d2 )

data = dlmread( 'crall_019_red.csv' )
d2 = data(2:end,:);
d2(:,2) = sqrt(d2(:,2));
d2(:,3) = sqrt(d2(:,3));
d2(:,4) = sqrt(d2(:,4));
dlmwrite( 'crall_019_red2.csv', d2 )