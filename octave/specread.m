function data = specread(filename,nscan)
% Reads scans from a SPEC file
% nscan = desired scan numbers, can be a vector, if not given reads them all
% Returns them as a structure data.y
% If only one scan read returns a matrix
% KH 30.03.98

if (nargin<2)
   nscan=-1;
end
% First open the file, if error occurs exit the function
[fid,msg] = fopen(filename,'r');
if fid == -1
   error(msg);
end
ifound=0;		% Counter for found scans
line=fgetl(fid);	% Reads the first line
while (ischar(line) && ~feof(fid))
    if length(line)>3;
        if (strcmp(line(1:2),'#S'))
            iscan=sscanf(line(3:length(line)),'%g');							% Scan number
            if (sum(iscan==nscan)>0 || nargin<2)								% Correct scan
                ifound=ifound+1;												% Add scan counter
                data(ifound).y=[];												% Initialize data variable
                line(1:2)='##';
                while (ischar(line) && ~strcmp(line(1:2),'#S'));				% Empty or scan end
                    line=fgetl(fid);											% Read new line
                    line=[line '  '];											% Trick to accept empty line
                    [numdata,n]=sscanf(line,'%g');
                    if n > 0;
                        numdata=numdata'; 
						if size(data(ifound).y,2)==0 || n==size(data(ifound).y,2)
							data(ifound).y=[data(ifound).y;numdata];	% Read numbers 
						end
                    end
                    if feof(fid)
                        break;
					end                    
                end
            else
                line=fgetl(fid);  
            end
        else
            line=fgetl(fid);
        end
    else
        line=fgetl(fid);
    end
end
fclose(fid);
if (nargin==2);
   if ifound~=length(nscan);
      warning('Not all scans found!');
   elseif length(nscan)==1;
      data = data(1).y;
   end
end
