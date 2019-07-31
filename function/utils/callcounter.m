function y = callcounter(A, x, mu)

    if ~isa(A, 'function_handle')
        error('A must be a function handle!\n')
        y = 0;
    end

    global calls;

    if nargin==2
        calls = calls + 1;
        y = A(x);
    else
        calls = calls + 1;
        y = A(x,mu);
    end
end

