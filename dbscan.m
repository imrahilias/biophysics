#!/bin/octave
## This reposetory contains a simple implementation of DBSCAN algorithm
## using GNU OCTAVE.  DBSCAN is a popular clustering algorithm. It
## creates clusters on a spatial data depending on two parameters:
## MinPoints: minimum number of points needed in its neighbourhood to
## consider it as a valid data(not noise).  EPS: A distance on which
## neighbourhood is calculated.
## For more info:
## https://github.com/devil1993/DBSCAN
## https://en.wikipedia.org/wiki/DBSCAN
## http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=1A6A7A85AF3F43BCBF66D847FEC8F8C5?doi=10.1.1.121.9220&rep=rep1&type=pdf
function [assignments,C] = dbscan(X,minpts,EPS)
  C = 0;
  assignments = zeros(size(X)(1),1);
  clustered = zeros(size(X)(1),1);
  for i=1: size(X)(1)
    if(clustered(i)==1)
      continue;
    endif
    clustered(i)=1;
    isneighbour = [];
    neighbourcount = 0;
    for j=1: size(X)(1)
      dist = sqrt(sum((X(i,:)-X(j,:)).^2));
      if(dist<EPS)
        neighbourcount++;
        isneighbour = [isneighbour j];
      endif
    endfor
    if(neighbourcount<minpts)
      continue;
    else
      C++;
      assignments(i) = C;
      for k=isneighbour
        if(clustered(k)==0)
          clustered(k) = 1;
          _isneighbour = [];
          _neighbourcount = 0;
          for j=1: size(X)(1)
            dist = sqrt(sum((X(k,:)-X(j,:)).^2));
            if(dist<EPS)
              _neighbourcount++;
              _isneighbour = [_isneighbour j];
            endif
          endfor
          if(_neighbourcount>=minpts)
            isneighbour = [isneighbour _isneighbour];
          endif
        endif
        assignments(k) = C;
      endfor
    endif
  endfor
