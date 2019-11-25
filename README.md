# disconnected-loop
2019-09-19:
update the xml input for displacements;
<Displacement>
	<link_max>1</link_max> \\default: 1
	<link_dirs>0 1 2 3 4 5 6 7</link_dirs> \\default: 0 1 2 3 4 5 6 7
</Displacement>


2019-11-20:
previous version is for LaMET-type operator.

update the support for second moment operator, covering all possible 2-links.

2019-11-23:

Do LaMET-type for link_max>2, do moments for link_max<3.

Output link pattern are named in different ways. For moment the numbers are {length}_{dir1}_{dir2}, for LaMET the numbers are {dir}_{length}.


2019-11-25:

Do first and second moments for all cases, do LaMET for link_max>2, save to different files.

