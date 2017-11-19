# See also
* [Computing external dynamic angle](https://gitlab.com/adammajewski/dynamic_external_angle)
* 

# git ( gitlab)

```
cd existing_folder
git init
git remote add origin git@gitlab.com:adammajewski/parameter_external_angle.git
git add .
git commit -m "Initial commit"
git push -u origin master
```


# technical note
GitLab uses:
* the Redcarpet Ruby library for [Markdown processing](https://gitlab.com/gitlab-org/gitlab-ce/blob/master/doc/user/markdown.md)
* KaTeX to render [math written with the LaTeX syntax](https://gitlab.com/gitlab-org/gitlab-ce/blob/master/doc/user/markdown.md), but [only subset](https://khan.github.io/KaTeX/function-support.html)

