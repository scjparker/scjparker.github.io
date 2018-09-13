#! /bin/sh

## edit if needed; where'e the [remote] for repository
REMOTE=origin
BRANCH=gh-pages

## Build: requires `bundler`
bundle exec jekyll build

## Add and commit build changes
git add -A
git commit -m "Build ${date}"

## git subtree deployment to gh-pages here
## basically, just pushes the `_site` directory to the specified <remote>/<branch>
git subtree push --prefix _site origin gh-pages
