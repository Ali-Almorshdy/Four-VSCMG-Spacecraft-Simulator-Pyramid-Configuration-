function s = MRPswitch(q)

% MRPswitch
%
%	S = MRPswitch(Q,s2) checks to see if norm(Q) is larger than s2.
%	If yes, then the MRP vector Q is mapped to its shadow set.
%

q2 = q'*q;
if(q2>=1)
	s = -q/q2;
else
	s = q;
end
