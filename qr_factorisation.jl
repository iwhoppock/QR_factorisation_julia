#Classical GS
function qrcgs(A)
  (m,n) = size(A)
  Q = zeros(m,n)
  R = zeros(n,n)
  for j=1:m
    aj = A[:,j]
    vj = aj
    for i = 1:j-1
      R[i,j] = dot(Q[:,i], aj)
      vj -= Q[:,i] * R[i,j]
    end
    R[j,j] = norm(vj)
    Q[:,j] = vj / R[j,j]
  end
    return Q,R
end;

#Modified GS
function qrmgs(A)
  (m,n) = size(A)
  Q = zeros(m,n)
  R = zeros(n,n)
  for j = 1:m
    aj = A[:,j]
    vj = aj
    for i = 1:j-1
      R[i,j] = dot(Q[:,i], vj)
      vj -= Q[:,i]*R[i,j]
    end
    R[j,j] = norm(vj)
    Q[:,j] = vj / R[j,j]
  end
    return Q,R
end;

#Householder
function qrhouse(A)
    (m,n) = size(A)
    Q = eye(m)
    R = float(A)
    for i = 1:n
        x = R[i:m,i]
        e = zeros(length(x))
        e[1] = 1
        vk = sign(x[1]) * norm(x) * e + x
        vk = vk / norm(vk)
        R[i:m,i:n] -= 2 * vk * vk' * R[i:m,i:n]
        Q[i:m,i:n] -= Q[i:m,i:n] * 2 * vk * vk'
    end
    return Q,R
end;

function backsub(R,b)
    n = length(b)
    x = zeros(n,1)
    for i = n:-1:1
        x[i] = b[i]
        for j = i+1:n
            x[i] = x[i] - R[i,j] * x[j]
        end
        x[i] = x[i] / R[i,i]
    end
    return x
end;
