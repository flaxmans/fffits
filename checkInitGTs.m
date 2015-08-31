function checkInitGTs()

% unit testing for setting up initial genotypes using initial allele counts

system('cut -f2 -d, InitialAlleleFreqs.csv | tail +2 > initCounts.txt');

initialCounts = importdata('initCounts.txt'); % aggregate allele count data (NOT from genotypes)

initialGTs = importdata('InitialGenotypes.csv');
igts = initialGTs.data;
igtCounts = sum(igts); % get total derived allele counts from genotypes

if ( 2*(numel(initialCounts)) ~= numel(igtCounts) )
    error('You suck; lengths off compared to what you expect');
end

countCheck = zeros(size(initialCounts));
n = numel(initialCounts);

for i = 1:n
    countCheck(i) = igtCounts((2*i)-1) + igtCounts(2*i); % sum the two haploid counts for one diploid count
end

if all(countCheck == initialCounts) && all(countCheck > 0)
    disp('InitialGenotypes.csv checks out against InitialAlleleFreqs.csv')
else
    error('InitialGenotypes.csv does NOT check out against InitialAlleleFreqs.csv')
end

system('rm initCounts.txt');
