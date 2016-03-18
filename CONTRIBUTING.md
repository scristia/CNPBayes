# Contributing

Follow these steps for adding a feature or fixing a bug. This will update the devel branch on Bioconductor

```
git checkout master
git merge branch-with-updates
git checkout devel
# look in git log to find the commit hashes
git cherry-pick oldest-new-commit..newest-new-commit 
git svn info # should just give repository information
git svn rebase # pulls down changes from svn repository, may be conflicts
git svn dcommit # "pushes" changes to svn
```
